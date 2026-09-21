############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
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
        residuals=peak.residuals.mx[peaks.to.keep, samples.to.keep],
        locations=peak.locations.gr[peaks.to.keep, ],
        sample.metadata=sample.metadata %>% filter(SampleID %in% samples.to.keep)
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
    # Run SVA 
    # full.fml %>% model.matrix(sample.metadata) %>% dim()
    # reduced.fml %>% model.matrix(sample.metadata) %>% dim()
    svs <- 
        tryCatch(
            {
                counts.matrix %>%
                as.matrix() %>%
                {.[rowSums(.) > 0, ]} %>%
                # svaseq(
                sva(
                    # full model matrix
                    full.fml %>% model.matrix(sample.metadata),
                    # covariate only model matrix
                    reduced.fml %>% model.matrix(sample.metadata),
                    # ...
                ) %>%
                {.$sv} %>%
                set_colnames(paste0('SV', 1:ncol(.))) %>% 
                as.data.frame() %>% 
                as_tibble()
            },
            error=function(cond){
                message('Estimating SVs failed with following error:')
                message(conditionMessage(cond))
            }
        )
}

run_sva_on_peak_residuals <- function(
    residuals.filepath,
    coords.filepath,
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    full.SV.model.vars=c('AD_CERAD_withDLB', 'Braak_3levels', 'CDR_3levels', 'age'),
    reduced.SV.model.vars=NULL,
    ...) {
    # subset + order data as defined
    peak.data <- 
        pick_samples_to_use(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # remove any rows with NAs in any model variable
    sample.metadata <- 
        sample.metadata %>%
        filter(SampleID %in% colnames(peak.data$residuals)) %>%
        filter(
            if_all(
                all_of(full.SV.model.vars),
                ~ !is.na(.)
            )
        )
    # remove any variables that are uniform across all samples
    for (model.var in full.SV.model.vars) {
        if (length(unique(sample.metadata[[model.var]])) == 1) {
            # sample.metadata <- sample.metadata %>% select(-c(!!sym(model.var)))
            full.SV.model.vars <- full.SV.model.vars[full.SV.model.vars != model.var]
            message(glue('removed {model.var} from SV model since it is singular'))
        }
    }
    # estimate SVs from peak residuals 
    run_sva(
        # counts.matrix=peak.data$residuals[1:1000, sample.metadata$SampleID],
        counts.matrix=peak.data$residuals[, sample.metadata$SampleID],
        sample.metadata=sample.metadata,
        full.model=full.SV.model.vars,
        reduced.model=reduced.SV.model.vars
    )
}

list_all_SVs <- function(){
    SVA_RESULTS_DIR %>% 
    parse_results_filelist(
        suffix='peak.residual.SVs.tsv',
        filename.column.name='filename'
    ) %>%
    select(-c(filename))
}

adjust_residuals_with_svs <- function(
############################################################
# Generate decorate clusters
############################################################
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
    cores,
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
    # tree.list %>% save_correlation_matrices(results_dir=results_dir)
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
        # scoreClusters(all.tree.list.clusters, BPPARAM=SnowParam(cores))
        scoreClusters(all.tree.list.clusters)
    # tree.list
    # all.tree.list.clusters
    scores.df <- 
        tree.scores %>%
        sapply(
            FUN=as_tibble,
            simplify=FALSE,
            USE.NAMES=TRUE
        ) %>% 
        bind_rows(.id='meanClusterSize')

    tree.scores$`10` %>% as_tibble()
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
        tree.list=tree.list,
        cluster.LEFs=,
        all.clusters=all.tree.list.clusters,
        filtered.clusters=
    )
}

run_decorate_pipeline <- function(
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

