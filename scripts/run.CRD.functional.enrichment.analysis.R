############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))

############################################################
# Set up input datasets
############################################################
# List all called CRD results
all.CRDs.df <- 
    list_all_results_files_in_set(set.name='merged.CRDs')
# load mapping of peaks to associated genes
# map all peal associated genes to their functions, so 
# I can get a list of all associated PeakIDs for each MSigDB term
all.MSigDB.terms.df <- 
    load_all_peak_metadata() %>% 
    filter(CPM.cutoff == '3CPM') %>%
    select(PeakID, gene.EnsemblID, gene.name) %>% 
    filter(!is.na(gene.name)) %>% 
    inner_join(
        load_all_msigdb_annotations() %>% unnest(term.genes),
        by=join_by(gene.name == term.genes)
    ) %>%
    group_by(annotation.souce, term.name) %>% 
    summarize(
        .groups='drop',
        across(
            .cols=c(gene.name, gene.EnsemblID, PeakID),
            .names='term.{.col}s',
            .fns=~ .x[!is.na(.x)] %>% unique() %>% list()
        )
    ) %>%
    group_by(annotation.souce) %>%
    mutate(
        # how many tests to correct for
        annotation.nterms=dplyr::n(),
        # size of the universe for fisher testing i.e. 
        # number of unique peaks mapped to any gene across all annotation terms
        annotation.total.peaks=term.PeakIDs %>% c() %>% unlist() %>% unique() %>% length(),
        # size of the universe for fisher testing i.e. 
        # number of unique genes across all annotation terms
        annotation.total.genes=term.gene.names %>% c() %>% unlist() %>% unique() %>% length()
    ) %>%
    ungroup()

all.MSigDB.terms.df %>% head(1) %>% t()
