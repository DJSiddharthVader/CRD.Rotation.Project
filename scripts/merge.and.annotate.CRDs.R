############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
all.sample.metadata <- load_sample_metadata()

############################################################
# List all generated CRD results
############################################################
plan(multisession, workers=TOTAL_CORES)
CRDs.df <- 
    list_all_results_files_in_set(set.name='CRDs') %>% 
    filter(n.SVs == NUM_SVS_TO_GENERATE) %>% 
    filter(SVs.included == NUM_SVS_TO_REMOVE) %>% 
    mutate(
        CRDs=
            future_pmap(
                .l=list(filepath),
                .f=~ read_tsv(.x, show_col_types=FALSE),
                .progress=TRUE
            ) 
    ) %>% 
    select(-c(filepath)) %>% 
        {.}; CRDs.df; CRDs.df$CRDs[[1]]
    CRDs.df %>% 
    unnest(CRDs) %>%
    count(residual.model, sample.strategy, is.decorate.filtered)

