############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(vegan)
    library(sva)
    library(GenomicRanges)
    library(decorate)
    library(limma)
    library(variancePartition)
})

############################################################
# Input PreProcessing
############################################################
select_sample_subset <- function(
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    ...){
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='AD.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='Control.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='No.Others'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='All.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='Make.Error'
    # decide which specific samples to select  based on inclusin criteria
    filter_cmd <- 
        if        (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'Control.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Control')
        } else if (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'AD.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'AD')
        } else if (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'No.Others') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) != 'Other')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'Control.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Braak_3to4')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'AD.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Braak5plus')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'No.Others') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) != 'BraakMax2')
        } else if (sample.strategy == 'All.Samples') {
            . %>% filter(!is.na(!!sym(AD.definition.column)))
        } else {
            stop(glue('sample inclusion criteria not implemented. you used {sample.strategy} to filter on the column {AD.definition.column}'))
        }
    # Now pull the Sample_IDs from the metadata for the selected samples
    sample.metadata %>%
    filter_cmd() %>% 
    pull(SampleID)
}

subset_peak_data <- function(
    residuals.filepath,
    coords.filepath,
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    ...){
    # load coordinates of all ATAC features
    peak.locations.gr <- 
        coords.filepath %>% 
        readRDS() %>% 
        {.$genes} %>% 
        as.data.frame() %>% 
        makeGRangesFromDataFrame()
        # select(seqnames, start, end)
    # load residual ATAC matrix
    peak.residuals.mx <- 
        residuals.filepath %>% 
        readRDS()
    # Only keep common peaks
    peaks.to.keep <- 
        intersect(
            names(peak.locations.gr),
            rownames(peak.residuals.mx)
        )
    # select which samples to include when annotating CRDs 
    samples.to.keep <- 
        select_sample_subset(
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # Now make sure rows are compatible
    # length(samples.to.keep); length(peaks.to.keep)
    # peak.locations.gr <- peak.locations.gr[peaks.to.keep, ]
    # peak.residuals.mx <- peak.residuals.mx[peaks.to.keep, samples.to.keep]
    list(
        sample.metadata=sample.metadata %>% filter(SampleID %in% samples.to.keep),
        residuals=peak.residuals.mx[peaks.to.keep, samples.to.keep],
        locations=peak.locations.gr[peaks.to.keep, ]
    )
}

############################################################
# Run SVA on peak residuals
############################################################
run_sva <- function(
    sample.metadata,
    counts.matrix,
    full.model.vars,
    reduced.model.vars=NULL,
    n.SVs='leek',
    z.score.counts=FALSE,
    ...) {
    # Run SVA to identify SVs and return SV matrix
    # First convert variable lists to formula objs
    reduced.fml <- 
        if (is.null(reduced.model.vars)) {
            formula('~ 1')
        } else {
            reduced.model.vars %>%
            paste(collapse='+') %>%
            sprintf('~%s', .)
            formula() 
        }
    # effects  + covaraites
    full.fml <- 
        full.model.vars %>%
        paste(collapse='+') %>%
        sprintf('~%s', .) %>% 
        formula() 
    # z-score transform features (rows) of the count matrix if specified
    counts <- 
        counts.matrix %>%
        as.matrix() %>%
        {.[rowSums(.) > 0, ]} %>%
        {
            if (z.score.counts) {
                t(scale(t(.), center = TRUE, scale = TRUE)) # compute z-score for each row (OCR) 
            } else {
                .
            }
        }
    # specify model matrices + number of SVs for SV estimation
    mod <- 
        full.fml %>% model.matrix(sample.metadata)
    mod0 <- 
        reduced.fml %>% model.matrix(sample.metadata)
    n.sv <- 
        case_when(
            n.SVs == 'leek' ~ num.sv(dat=counts, mod=mod, method='leek'),
            n.SVs == 'be'   ~ num.sv(dat=counts, mod=mod, method='be'),
            TRUE            ~ as.integer(n.SVs)
        )
    message(full.fml)
    message(reduced.fml)
    # message(c(n.SVs, n.sv))
    # Run SVA to estimate SVs as specified from the count matrix
    sva(
        dat=counts,
        mod=mod,
        mod0=mod0,
        n.sv=n.sv
    ) %>%
    {.$sv} %>%
    set_colnames(paste0('SV', 1:ncol(.))) %>% 
    as.data.frame() %>% 
    as_tibble() %>%
    add_column(SampleID=colnames(counts.matrix))
}

run_sva_on_peak_residuals <- function(
    residuals.filepath,
    coords.filepath,
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    full.SV.model.vars=c('AD_CERAD_withDLB', 'Braak_3levels', 'CDR_3levels', 'age'),
    reduced.SV.model.vars=NULL,
    n.SVs='leek',
    z.score.counts=FALSE,
    ...) {
    # subset + order data as defined
    peak.data <- 
        subset_peak_data(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # remove any rows with NAs in any model variable
    all.vars <- c(full.SV.model.vars, reduced.SV.model.vars)
    sample.metadata <- 
        sample.metadata %>%
        filter(SampleID %in% colnames(peak.data$residuals)) %>%
        mutate(across(where(is.factor) & all_of(all.vars), droplevels)) %>%
        # Need to remove rows with NAs to make matrix solvable
        filter(!if_any(all_of(all.vars), is.na))
    # remove any variables that are uniform across all samples
    for (model.var in c(full.SV.model.vars, reduced.SV.model.vars)) {
        if (length(unique(sample.metadata[[model.var]])) == 1) {
            # sample.metadata <- sample.metadata %>% select(-c(!!sym(model.var)))
            full.SV.model.vars <- full.SV.model.vars[full.SV.model.vars != model.var]
            message(glue('removed {model.var} from SV model since it is uniform across samples'))
        }
    }
    # estimate SVs from peak residuals 
    run_sva(
        counts.matrix=peak.data$residuals[, sample.metadata$SampleID],
        sample.metadata=sample.metadata,
        full.model=full.SV.model.vars,
        reduced.model=reduced.SV.model.vars,
        n.SVs=n.SVs,
        z.score.counts=z.score.counts,
    )
}

residualize_out_specified_SVs <- function(
    peak.residuals.mx,
    SVs.df,
    SVs.included=NULL){
    # peak.residuals.mx=peak.data$residuals; SVs.df=SVs.df; SVs.included=NULL
    if (is.null(SVs.included)) { SVs.included <- c(colnames(SVs.df)) }
    model.formula <- 
        SVs.df %>% 
        colnames() %>% 
        # {.[1:(SVs.included)]} %>% 
        {.[. == SVs.included]} %>% 
        paste(collapse="+") %>% 
        sprintf('~ %s', .) %>% 
        formula()
    common.samples <- 
        intersect(rownames(SVs.df), colnames(peak.residuals.mx))
    # Regress out SVs
    # peak.residuals.mx[, common.samples, drop=FALSE] %>% 
    peak.residuals.mx[, common.samples] %>% 
    lmFit(
        model.matrix(
            model.formula,
            SVs.df[common.samples, , drop=FALSE]
        )
    ) %>%
    residuals(peak.residuals.mx[, common.samples])
    # residuals(peak.residuals.mx[, common.samples, drop=FALSE])
}

residualize_SVs <- function(
    SVs.included,
    SVs.df,
    residuals.mx){
    if (SVs.included == 0) {
        message('Not regressing out any SVs before running decorate')
        residuals.mx
    } else if (SVs.included > 0) { 
        SV.variables <- 
            SVs.df %>% colnames() %>% grep('^SV', ., value=TRUE)
        SVs.included <- 
            min(SVs.included, length(SV.variables))
        message(glue('Regressing out {SVs.included} SVs before running decorate'))
        residualize_out_specified_SVs(
            peak.residuals.mx=residuals.mx,
            SVs.df=SVs.df,
            SVs.included=SV.variables[1:SVs.included]
        )
    } else {
        stop(glue('Invalid arg to SVs.included: {SVs.included}'))
    }
}

generate_elbow_data <- function(
    SVs.filepath,
    all.sample.metadata,
    BPPARAM=SerialParam(),
    ...){
    # paste0(c('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; ')), collapse='; ')
    all.SVs <- 
        SVs.filepath %>% 
        read_tsv(show_col_types=FALSE, n_max=1) %>% 
        colnames()
    tibble(SVs.included=0:length(all.SVs)) %>%
    mutate(
        LEFs=
            future_pmap(
                .l=.,
                .f=run_decorate_pipeline,
                ...,
                SVs.filepath=SVs.filepath,
                sample.metadata=all.sample.metadata,
                LEFs.only=TRUE,
                p=progressor(label='SVs included', steps=length(all.SVs)),
                .progress=TRUE
            )
    ) %>%
    unnest(LEFs)
}

run_varaincePartition <- function(
    peak.residuals.mx,
    sample.metadata,
    model.formula,
    SVs.included=0,
    BPPARAM=SerialParam(),
    p=NULL,
    ...){
    # peak.residuals.mx=peak.data$residuals
    common.samples <- 
        colnames(peak.residuals.mx) %>%
        intersect(rownames(sample.metadata))
    # residualize SVs out of input matrix before running variancePartition()
    peak.residuals.mx <-
        residualize_SVs(
            SVs.included=SVs.included,
            SVs.df=sample.metadata[common.samples, ],
            residuals.mx=peak.residuals.mx[, common.samples]
        )
    # print formula 
    message(model.formula)
    # estiamte how much each formula variable contributes to variance in each 
    # peak's abundance in the residualized matrix
    vps.df <- 
        fitExtractVarPartModel(
            form=model.formula,
            exprObj=peak.residuals.mx[, common.samples],
            data=sample.metadata[common.samples, ],
            ...
        ) %>% 
        rownames_to_column('PeakID') %>% 
        as_tibble()
    if (!is.null(p)) { p() }
    return(vps.df)
}

generate_SV_varpar_results <- function(
    all.sample.metadata,
    residuals.filepath,
    coords.filepath,
    AD.definition.column,
    sample.strategy,
    SVs.filepath,
    model.variables,
    ...){
    # paste0(c('row.index=1; model.variables=RELEVANT_METADATA_COLUMNS; BPPARAM=SnowParam(TOTAL_CORES * (4 / 5));', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; ')), collapse='; ')
    # load + subset data according to specified params
    peak.data <- 
        subset_peak_data(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=all.sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # join SVs with sample metadata for variancePartition
    sample.metadata <- 
        # Load SVs to residualize out
        SVs.filepath %>%
        read_tsv(show_col_types=FALSE) %>%
        as.data.frame() %>%
        # match SVs with sample metadata + phenotypes
        merge(as.data.frame(peak.data$sample.metadata), by='SampleID') %>%
        column_to_rownames('SampleID')
    # run variancePartition with the formula ~ model.variables + SV1 + ... + SVi
    # build model formula automatically from provided SVs
    SV.variables <- 
        sample.metadata %>% colnames() %>% grep('^SV', ., value=TRUE)
    model.formula <- 
        # metadata to include in the formula for variancePartition
        pmap_lgl(
            .l=list(model.variables),
            # automatically remove any variables that are uniform across all samples 
            .f=~ length(unique(sample.metadata[[.x]])) != 1
        ) %>% 
        {model.variables[.]} %>% 
        c(SV.variables) %>% 
        paste0(collapse='+') %>% 
        sprintf('~ %s', .) %>% 
        formula()
    # want to compare how variances explained for model.variable changes across 
    # different number of SVs incorporated, also a kind of elbow analysis
    # tmp[row.index, ] %>% cross_join(tibble(SVs.included=0:length(SV.variables))) %>% 
    tibble(SVs.included=0:length(SV.variables)) %>%
    mutate(
        variancePartitionResults=
            future_pmap(
                .l=.,
                .f=run_varaincePartition,
                peak.residuals.mx=peak.data$residuals,
                ..., 
                model.formula=model.formula,
                sample.metadata=sample.metadata,
                hideErrorsInBackend=TRUE, showWarnings=TRUE,
                p=progressor(label='SVs included', steps=length(SV.variables)),
                .progress=TRUE
            ) 
    ) %>% 
    add_column(
        total.SVs=length(SV.variables),
        model.variables=paste0(model.variables, collapse=', ')
    ) %>% 
    unnest(variancePartitionResults)
}

############################################################
# Generate decorate clusters
############################################################
generate_LEFs_only_with_decorate <- function(
    peak.residuals.mx,
    peak.locations,
    adjacentCount,
    method.corr,
    clusterMethod,
    meanClusterSize,
    BPPARAM=SerialParam(),
    ...){
    # Evaluate hierarchical clustering
    # adjacentCount is the number of adjacent peaks considered in correlation
    # use Spearman correlation to reduce the effects of outliers
    # peak.residuals.mx=peak.data$residuals; peak.locations=peak.data$locations;
    tree.list <- 
        peak.residuals.mx %>% 
        runOrderedClusteringGenome( 
            peak.locations,
            adjacentCount=adjacentCount,
            method.corr=method.corr
        )
    # Choose cutoffs and return clusters using multiple values for meanClusterSize 
    # Clusters corresponding to each parameter value are returned and then processed downstream
    # By using multiple parameter values, epigenetic features are included in clusters 
    all.tree.list.clusters <- 
    # at different resolutions
        tree.list %>% 
        createClusters(
            method=clusterMethod,
            meanClusterSize=meanClusterSize
        )
    # return clusters + LEFs across meanClusterSizes as tibble
    tree.list %>% 
    scoreClusters(
        all.tree.list.clusters,
        BPPARAM=BPPARAM
    ) %>% 
    sapply(
        FUN=as_tibble,
        simplify=FALSE,
        USE.NAMES=TRUE
    ) %>% 
    bind_rows(.id='meanClusterSize')
}

generate_peak_cluster_with_decorate <- function(
    peak.residuals.mx,
    peak.locations,
    adjacentCount,
    method.corr,
    clusterMethod,
    meanClusterSize,
    jaccardCutoff,
    filterMetric,
    filterMetricCutoff,
    BPPARAM=SerialParam(),
    ...){
    # Evaluate hierarchical clustering
    # adjacentCount is the number of adjacent peaks considered in correlation
    # use Spearman correlation to reduce the effects of outliers
    # peak.residuals.mx=peak.data$residuals; peak.locations=peak.data$locations;
    tree.list <- 
        peak.residuals.mx %>% 
        runOrderedClusteringGenome( 
            peak.locations,
            adjacentCount=adjacentCount,
            method.corr=method.corr
        )
    # Choose cutoffs and return clusters using multiple values for meanClusterSize 
    # Clusters corresponding to each parameter value are returned and then processed downstream
    # By using multiple parameter values, epigenetic features are included in clusters 
    # at different resolutions
    all.tree.list.clusters <- 
        tree.list %>% 
        createClusters(
            method=clusterMethod,
            meanClusterSize=meanClusterSize
        )
    tree.scores <- 
        tree.list %>% 
        scoreClusters(
            all.tree.list.clusters,
            BPPARAM=BPPARAM
        )
    filtered.tree.list.clusters <- 
        tree.scores %>% 
        retainClusters( 
            metric=filterMetric,
            cutoff=filterMetricCutoff
        ) %>% 
        # get retained clusters
        filterClusters(all.tree.list.clusters, .) %>% 
        # collapse redundant clusters
        collapseClusters(
            peak.locations,
            jaccardCutoff=jaccardCutoff
        )
    # return all data
    list(
        cluster.LEFs=
            tree.scores %>%
            sapply(
                FUN=as_tibble,
                simplify=FALSE,
                USE.NAMES=TRUE
            ) %>% 
            bind_rows(.id='meanClusterSize'),
        tree.list=tree.list,
        all.clusters=all.tree.list.clusters,
        filtered.clusters=filtered.tree.list.clusters
    )
}

run_decorate_pipeline <- function(
    residuals.filepath,
    coords.filepath,
    SVs.filepath,
    sample.metadata,
    AD.definition.column='AD_CERAD_withDLB',
    sample.strategy='All',
    SVs.included=Inf,
    LEFs.only=FALSE,
    p=NULL,
    ...) {
    # subset + order data as defined
    peak.data <- 
        subset_peak_data(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # Load SVs to residualize out 
    if (SVs.included > 0) {
        SVs.df <- 
            SVs.filepath %>%
            read_tsv(show_col_types=FALSE) %>%
            as.data.frame() %>%
            column_to_rownames('SampleID')
        peak.residuals.mx <-
            residualize_SVs(
                SVs.df=SVs.df,
                residuals.mx=peak.data$residuals,
                SVs.included=SVs.included
            )
    } else  if (SVs.included == 0) {
        peak.residuals.mx <- 
            peak.data$residuals
    } else {
        stop(glue('SVs.included must be an integer, passed: {SVs.included}'))
    }
    # generate peak clusters, either just the cluster LEF scores or the full blob with peak-clusters
    if (LEFs.only) {
        results.obj <- 
            generate_LEFs_only_with_decorate(
                peak.residuals.mx=peak.residuals.mx,
                peak.locations=peak.data$locations,
                ...
            )
    } else {
        results.obj <- 
            generate_peak_cluster_with_decorate(
                peak.residuals.mx=peak.residuals.mx,
                peak.locations=peak.data$locations,
                ...
            )
    }
    if (!is.null(p)) { p() }
    return(results.obj)
}

############################################################
# Misc posterity code
############################################################
sva_docs_example_correlations <- function() {
    library(bladderbatch)
    data(bladderdata)
    pheno = pData(bladderEset)
    edata = exprs(bladderEset)
    mod = model.matrix(~as.factor(cancer), data=pheno)
    mod0 = model.matrix(~1,data=pheno)
    n.sv = num.sv(edata,mod,method="be"); n.sv
    svobj = sva(edata,mod,mod0,n.sv=n.sv)
    colnames(svobj$sv) <- paste0('SV', 1:n.sv)
    canCorPairs(
        paste(c('~ cancer', colnames(svobj$sv)), collapse='+'),
        data=bind_cols(svobj$sv, pheno)
    )
    n.sv = num.sv(edata,mod,method="leek"); n.sv
    svobj = sva(edata,mod,mod0,n.sv=n.sv)
    colnames(svobj$sv) <- paste0('SV', 1:n.sv)
    canCorPairs(
        paste(c('~ cancer', colnames(svobj$sv)), collapse='+'),
        data=bind_cols(svobj$sv, pheno)
    )
    # my function shich should be the same
    run_sva(
        sample.metadata=as_tibble(pheno ),
        counts.matrix=edata,
        full.model.vars=c('cancer'),
        reduced.model.vars=NULL,
    ) %>% 
    as.data.frame() %>% 
    column_to_rownames('SampleID') %>% 
    {
        canCorPairs(
            paste(c('~ cancer', colnames(.)), collapse='+'),
            data=bind_cols(., pheno)
        )
    }
}

