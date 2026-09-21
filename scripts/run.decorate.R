############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
sample.metadata <- load_sample_metadata()

############################################################
# Generate all input dataset + hyper-param combinations to test
############################################################
all.input.and.parameter.combinations.df <- 
    # First list residual ATAC-seq count matrices that can be used to annotate CRDs
    list_all_ATAC_residual_sets() %>% 
    # get data specifying peak coordinates for each set of residuals
    left_join(
        list_all_ATAC_coordinates_sets(),
        by=join_by(CPM.cutoff)
    ) %>% 
    filter(CPM.cutoff == '3CPM') %>% 
    filter(residual.model == 'KeepDxAge') %>% 
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
# Compute Surrogate Variables (SVs) and produce corrected matrices
############################################################
all.input.and.parameter.combinations.df %>% 
    mutate(
        results_dir=
            file.path(
                SVA_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('sample.strategy_{sample.strategy}')
            ),
        results_file=file.path(results_dir, 'peak.residual.SVs.tsv')
    ) %>% 
        # {.} -> tmp; tmp
        # tmp %>% 
    pmap(
        .f=check_cached_results,
        results_fnc=run_sva_on_peak_residuals,
        # .f=run_sva_on_peak_residuals,
        sample.metadata=sample.metadata,
        full.SV.model.vars=c('AD_CERAD_withDLB', 'Braak_3levels', 'CDR_3levels', 'age'),
        reduced.SV.model.vars=NULL,
    )
# list_all_SVs()

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
