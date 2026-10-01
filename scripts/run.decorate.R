############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
all.sample.metadata <- load_sample_metadata()

############################################################
# Generate all input dataset + hyper-param combinations to test
############################################################
all.input.and.parameter.combinations.df <- 
    list_all_ATAC_datasets() %>% 
    filter(CPM.cutoff == '3CPM') %>% 
    # which sets of samples to use to call CRDs
    cross_join(RESIDUAL_MATRIX_SAMPLE_PARAMS_DF) %>% 
    # map files with estimated SVs to the corresponding residual matrices
    inner_join(
        list_all_results_files_in_set(set.name='SVs'),
        by=
            join_by(
                CPM.cutoff,
                residual.model,
                sample.strategy,
                AD.definition.column
            )
    ) %>% 
    filter(z.score.counts, n.SVs == NUM_SVS_TO_GENERATE) %>% 
    cross_join(tibble(SVs.included=c(0,5))) %>%
    cross_join(DECORATE_HYPER_PARAMS_DF)

############################################################
# Run decorate for each parameter combination
############################################################
# For every residual matrix + sample strategy + hyper-param option 
# use decorate to annotate a set of CRDs
# all.input.and.parameter.combinations.df
N_PROCESSES <- TOTAL_CORES * (3 / 5)
N_CORES_PER_PROCESS <- floor((TOTAL_CORES - N_PROCESSES) / N_PROCESSES)
plan(multisession, workers=N_PROCESSES)
all.input.and.parameter.combinations.df %>% 
    mutate(
        results_dir=
            file.path(
                CRD_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('sample.strategy_{sample.strategy}'),
                glue('z.score.counts_{z.score.counts}'),
                glue('SVs.included_{SVs.included}'),
                glue('adjacentCount_{adjacentCount}'),
                glue('method.corr_{method.corr}'),
                glue('clusterMethod_{clusterMethod}'),
                glue('filterMetric_{filterMetric}'),
                glue('filterMetricCutoff_{filterMetricCutoff}'),
                glue('jaccardCutoff_{jaccardCutoff}')
            ),
        results_file=file.path(results_dir, 'decorate.blob.rds')
    ) %>% 
        {.} -> tmp; tmp
        # tmp %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        # force_redo=TRUE,
        # return_data=FALSE,
        # silence=TRUE,
        results_fnc=run_decorate_pipeline,
        sample.metadata=all.sample.metadata,
        BPPARAM=SnowParam(N_CORES_PER_PROCESS),
        p=progressor(label='decorate runs', steps=nrow(all.input.and.parameter.combinations.df)),
        .progress=TRUE
    ) 

