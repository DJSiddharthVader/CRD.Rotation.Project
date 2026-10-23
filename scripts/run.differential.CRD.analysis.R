############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))

############################################################
# Set up input datasets
############################################################
# load sample metadata
all.sample.metadata <- 
    load_sample_metadata()
# ATAC data
all.ATAC.and.SV.datasets.df <- 
    list_all_ATAC_datasets() %>% 
    # which sets of samples to use to call CRDs
    cross_join(RESIDUAL_MATRIX_SAMPLE_PARAMS_DF) %>% 
    # map files with estimated SVs to the corresponding residual matrices
    inner_join(
        list_all_results_files_in_set(set.name='SVs'),
        by=
            join_by(
                CPM.cutoff,
                residual.model,
                AD.definition.column,
                sample.strategy
            )
    ) %>% 
    filter(sample.strategy %in% c('All.Samples')) %>%
    select(-c(AD.definition.column, sample.strategy))
# Match CRDs to test with corresponding ATAC matrix + SVs to use as input for differential testing
all.input.paramater.combinations.df <- 
    # CRDs to test for differential signal
    list_all_results_files_in_set(set.name='merged.CRDs') %>% 
    # Match to relevant set of peaks 
    inner_join(
        all.ATAC.and.SV.datasets.df,
        by=
            join_by(
                CPM.cutoff,
                residual.model,
                z.score.counts,
                n.SVs
            )
    ) %>% 
    filter(CPM.cutoff == '3CPM') %>%
    filter(n.SVs == NUM_SVS_TO_GENERATE) %>%
    filter(z.score.counts) %>%
    filter(SVs.included == NUM_SVS_TO_REMOVE) %>% 
    filter(is.decorate.filtered == 'filtered') %>% filter(meanClusterSize == 'merged') %>% 
    # differential contrasts to test with decorate 
    # cross_join(DIFFERENTIAL_DECORATE_HYPER_PARAMS_DF %>% head(1))
    cross_join(DIFFERENTIAL_DECORATE_HYPER_PARAMS_DF)

############################################################
# Generate Differential CRD results
############################################################
plan(multisession, workers=TOTAL_CORES)
all.input.paramater.combinations.df %>% 
    mutate(
        results_dir=
            file.path(
                DIFFERENTIAL_CRD_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('CRD.sample.strategy_{sample.strategy}'),
                glue('z.score.counts_{z.score.counts}'),
                glue('n.SVs_{n.SVs}'),
                glue('SVs.included_{SVs.included}'),
                glue('adjacentCount_{adjacentCount}'),
                glue('method.corr_{method.corr}'),
                glue('clusterMethod_{clusterMethod}'),
                glue('filterMetric_{filterMetric}'),
                glue('filterMetricCutoff_{filterMetricCutoff}'),
                glue('jaccardCutoff_{jaccardCutoff}'),
                glue('meanClusterSize_{meanClusterSize}'),
                glue('is.decorate.filtered_{is.decorate.filtered}'),
                glue('contrast_{contrast}')
            ),
        results_file=file.path(results_dir, 'differential.CRD.results.tsv')
    ) %>% 
    # future_pmap(
    pmap(
        .l=.,
        .f=check_cached_results,
        # force_redo=TRUE,
        return_data=FALSE,
        results_fnc=compute_all_differential_CRD_results,
        all.sample.metadata=all.sample.metadata,
        test.methods=c('Delaneau', 'deltaSLE', 'Box.permute', 'Box'),
        .progress=TRUE
    ) 

