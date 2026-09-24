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
        # list_all_SVs() %>% head(1) %>% t()
        # list_all_SVs()
all.input.and.parameter.combinations.df <- 
    list_all_ATAC_datasets() %>% 
    filter(CPM.cutoff == '3CPM') %>% 
    # which sets of samples to use to call CRDs
    cross_join(RESIDUAL_MATRIX_SAMPLE_PARAMS_DF) %>% 
    # map files with estimated SVs to the corresponding residual matrices
    left_join(
        list_all_SVs(),
        by=join_by(CPM.cutoff, residual.model, AD.definition.column, sample.strategy)
    ) %>% 
    # Define all relevant arguments/analysis decisions/parameters to decorate in the results filepath
    # explicity, this also allows for easier parsing of the results, associating each results
    # file with this metadata automatically in R
    cross_join(DECORATE_HYPER_PARAMS_DF) %>%
    # cross_join(tibble(SVs.included=c('All', 'None', 'uncorrelated')))
    # cross_join(tibble(SVs.included=c('All', 'None', 'uncorrelated', 5)))
    cross_join(tibble(SVs.included=c('All', 'None')))

############################################################
# Run decorate for each parameter combination
############################################################
# For every residual matrix + sample strategy + hyper-param option 
# use decorate to annotate a set of CRDs
# all.input.and.parameter.combinations.df
all.input.and.parameter.combinations.df %>% 
    # Define all relevant arguments/analysis decisions/parameters to decorate in the results filepath
    # explicity, this also allows for easier parsing of the results, associating each results
    # file with this metadata automatically in R
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
            )
    ) %>% 
        {.} -> tmp; tmp
        # tmp %>% 
    pmap(
        .l=.,
        .f=generate_decorate_peak_clusters,
        # .f=check_cached_results,
        # results_fnc=generate_decorate_peak_clusters,
        sample.metadata=sample.metadata
    ) 

############################################################
# Load decorate CRD annotations
############################################################
decorate.results.df <- 
    list_all_decorate_peak_clusters() %>% 
    filter(cluster.scope == 'filtered.clusters')
# decorate.results.df %>% 
#     count(CPM.cutoff, residual.model, sample.strategy, meanClusterSize, cluster.scope)
# decorate.results.df$filepath[[1]] %>% read_tsv() %>% count(chr, CRDID) %>% count(chr)
