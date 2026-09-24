############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# handlers(global=TRUE) # for progress bars
# load patient metadata
all.sample.metadata <- 
    load_sample_metadata()
# List all input residual matrices 
peak.matrices.df <- 
    list_all_ATAC_datasets() %>% 
    filter(CPM.cutoff == '3CPM') %>% 
    # filter(residual.model == 'KeepDxAge') %>% 
    # hyper-params to select sample subsets for estimating SVs from 
    cross_join(RESIDUAL_MATRIX_SAMPLE_PARAMS_DF) %>%
    cross_join(
        expand_grid(
            n.SVs=c('leek', 1:NUM_SVS_TO_GENERATE),
            z.score.counts=c(TRUE, FALSE)
        )
    )

############################################################
# Compute Surrogate Variables (SVs) and produce corrected matrices
# generate with different numbers?
############################################################
plan(multisession, workers=TOTAL_CORES)
peak.matrices.df %>% 
    mutate(
        results_dir=
            file.path(
                SVA_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('AD.definition.column_{AD.definition.column}'),
                glue('sample.strategy_{sample.strategy}'),
                glue('n.SVs_{n.SVs}'),
                glue('z.score.counts_{z.score.counts}')
            ),
        results_file=file.path(results_dir, 'peak.residual.SVs.tsv')
    ) %>% 
    future_pmap(
        .f=check_cached_results,
        force_redo=TRUE,
        return_data=FALSE,
        results_fnc=run_sva_on_peak_residuals,
        sample.metadata=all.sample.metadata,
        full.SV.model.vars=ALL_PHENOTYPE_COLUMNS,
        reduced.SV.model.vars=NULL,
        .progress=TRUE
    )

############################################################
# Compute LEFs from residuals with varying numbers of SVs regressed out
peak.matrices.df %>% 
############################################################
# take all input residual matrices
# map all SV sets to residual matrices they can be applied to + decorate hyper-params
peak.and.SV.combos.df <- 
    peak.matrices.df %>% 
    # match them to all sets of SVs produced
    inner_join(
        list_all_SVs(),
        by=
            join_by(
                CPM.cutoff,
                residual.model,
                sample.strategy,
                AD.definition.column,
                n.SVs,
                z.score.counts
            )
    ) %>%
    # SVs are correlated no matter how many you estimate i.e. A=SVs[1..N] and B=SVs[1..M]
    # estimated from two different calls to run_sva() on the same input matrix
    # so if I=min(N,M) then cor(A[1..I], B[1..I] ~= 1 so 
    # we can just use the single largest set of SV estimates and vary how many SVs we residualize out ({1..N},{1..M})
    # to titrate the effects of SVs (instead of also measuring across SV sets ({A,B})
    filter(n.SVs == NUM_SVS_TO_GENERATE) %>% 
    # match against all combinations of hyper-params for decorate
    cross_join(DECORATE_HYPER_PARAMS_DF)
# Now iteratively residualize out SVs in integer order for each residual matrix + SVs then
# run decorate and save CRD LEFs to a file across various meanClusterSizes
plan(multisession, workers=TOTAL_CORES)
peak.and.SV.combos.df %>% 
    # define output filepath for cluster LEFs
    left_join(list_all_SVs()) %>%
    # for each SV i, include up to SVs 1:i and regress out
    # then run decorate to get CRD peak clusters + scores (LEFs)
    mutate(
        results_dir=
            file.path(
                ELBOW_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('AD.definition.column_{AD.definition.column}'),
                glue('sample.strategy_{sample.strategy}'),
                glue('n.SVs_{n.SVs}'),
                glue('z.score.counts_{z.score.counts}'),
                glue('adjacentCount_{adjacentCount}'),
                glue('method.corr_{method.corr}'),
                glue('clusterMethod_{clusterMethod}'),
                glue('filterMetric_{filterMetric}'),
                glue('filterMetricCutoff_{filterMetricCutoff}'),
                glue('jaccardCutoff_{jaccardCutoff}')
            ),
        # results_file=file.path(results_dir, 'elbow.cluster.data.tsv')
        results_file=file.path(results_dir, 'elbow.cluster.data.tsv')
    ) %>% 
    # for each SV i, include up to SVs 1:i and regress out -> call decorate clusters -> save LEFs
    # pmap(
    future_pmap(
        .l=.,
        .f=check_cached_results,
        return_data=FALSE,
        # force_redo=TRUE,
        silence=TRUE,
        results_fnc=generate_elbow_data,
        all.sample.metadata=all.sample.metadata,
        cores=1,
        .progress=TRUE

    )

############################################################
    )
