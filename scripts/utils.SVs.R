############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(sva)
    library(variancePartition)
})

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

############################################################
# Residualize out SVs from peaks
############################################################
residualize_out_specified_SVs <- function(
    peak.residuals.mx,
    SVs.df,
    SVs.included=NULL){
    if (is.null(SVs.included)) { SVs.included <- c(colnames(SVs.df)) }
    # Formula with only included SVs
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

############################################################
# Aanalyze SV-residualized matrices
############################################################
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
    residuals.filepath,
    all.sample.metadata,
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

regress_SVs_and_compute_stats <- function(
    peak.residuals.mx,
    sample.metadata,
    SVs.included,
    ...) { 
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
    # residual matrix to peak stats
    peak.residuals.mx %>% 
    as.data.frame() %>% 
    rownames_to_column('PeakID') %>%
    as_tibble() %>% 
    pivot_longer(
        -c(PeakID),
        names_to='SampleID',
        values_to='peak.residuals'
    ) %>% 
    group_by(PeakID) %>%
    summarize(
        .groups='drop',
        peak.variances=var(peak.residuals)
    )
}

compute_SV_regressed_peak_stats <- function(
    all.sample.metadata,
    residuals.filepath,
    coords.filepath,
    AD.definition.column,
    sample.strategy,
    SVs.filepath,
    ...){
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
    # want to compute simple summary stats for each peak for each 
    # different number of SVs regressed, also a kind of elbow analysis
    # tmp[row.index, ] %>% cross_join(tibble(SVs.included=0:length(SV.variables))) %>% 
    tibble(SVs.included=0:length(SV.variables)) %>%
    mutate(
        peak.stats=
            future_pmap(
                .l=.,
                .f=regress_SVs_and_compute_stats,
                peak.residuals.mx=peak.data$residuals,
                sample.metadata=sample.metadata,
                p=progressor(label='SVs included', steps=length(SV.variables)),
                .progress=TRUE
            ) 
    ) %>% 
    add_column(total.SVs=length(SV.variables)) %>% 
    unnest(peak.stats)
}

