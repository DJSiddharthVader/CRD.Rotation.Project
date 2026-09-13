############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'locations.R'))
source(here('scripts', 'utils.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
sample.metadata <- load_sample_metadata()


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

all.residual.sets.df %>%
    pmap(
        .l=.,
        .f=run_decorate,
        resids.matrix,
        resids.granges,
        adjacentCount=500,
        method.corr="spearman",
        clusterMethod="meanClusterSize",
        meanClusterSize=c(20, 30, 40, 50),
        filterMetric="LEF",
        filterMetricCutoff=0.15,
        jaccardCutoff=0.9
    ) 
