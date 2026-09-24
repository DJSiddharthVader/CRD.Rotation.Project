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
    left_join(list_all_SVs()) %>%
    # for each SV i, include up to SVs 1:i and regress out
    # then run decorate to get CRD peak clusters + scores (LEFs)
    mutate(
        results_dir=
            file.path(
                SVA_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('AD.definition.column_{AD.definition.column}'),
                glue('sample.strategy_{sample.strategy}')
            ),
        results_file=file.path(results_dir, 'elbow.plot.data.tsv')
    ) %>% 
    cross_join(DECORATE_HYPER_PARAMS_DF) %>% 
    add_column(cores=TOTAL_CORES) %>% 
        # {.} -> tmp; tmp
        # tmp %>% head(1) %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        # force_redo=TRUE,
        return_data=FALSE,
        results_fnc=generate_elbow_data,
        sample.metadata=all.sample.metadata,
        .progress=TRUE
    )
