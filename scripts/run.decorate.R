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
    # decorate hyper-parameter values to test
    cross_join(
        tibble(
            # Different inputs to generate CRD annotations from
            sample.strategy=
                c(
                    'Control.Samples',   # only annotated CRDs using Control samples
                    'AD.Samples',        # only annotate CRDs using AD samples
                    'All.Samples'        # annotate CRDs using All samples
                ),
            # which metadata column we use to define sample conditions
            AD.definition.column=
                c(
                    # 'Braak_3levels',   # Braak_5plus, Braak_3to4, Braak_Max2
                    'AD_CERAD_withDLB' # AD, Control, Other
                )
        )
    ) %>% 
    cross_join(
        tibble(
            adjacentCount=c(500),
            method.corr=c("spearman"),
            clusterMethod=c("meanClusterSize"),
            # meanClusterSize=list(list(20, 30, 40, 50)),
            meanClusterSize=list(list(10, 25, 50, 80, 100)),
            filterMetric=c("LEF"),
            filterMetricCutoff=c(0.15),
            jaccardCutoff=c(0.9)
        )
    )

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
