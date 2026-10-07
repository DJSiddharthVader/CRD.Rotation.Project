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

