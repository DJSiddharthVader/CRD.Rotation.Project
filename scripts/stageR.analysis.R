############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
library(stageR)

############################################################
# Set up input datasets
############################################################
# load sample metadata
all.sample.metadata <- load_sample_metadata()
# Load Differential OCR data
all.peak.differential.results.df <- 
    load_all_differential_OCR_results() %>%
    # load_all_differential_OCR_results(force_redo=TRUE) %>%
    filter(
        contrast %in% 
            c(
                'AD_CERAD_withDLB',
                'CDR_3levels',
                'Braak_3levels',
                'ApoE4count_2levels'
            )
    ) %>% 
    select(
        # CPM.cutoff, sample.strategy, residual.model,
        CPM.cutoff, residual.model,
        # chr, peak.start, peak.end, peak.length, 
        contrast, contrast.levels,
        PeakID, logFC, p.value, p.adjust
    )
all.merged.CRDs.df <- 
    list_all_results_files_in_set(set.name='merged.CRDs') %>% 
    filter(merging.level != 'separate') %>% 
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                .f=combine_CRD_Locations_and_PeakLists,
                .progress=TRUE
            )
    ) %>% 
    unnest(CRDs) %>% 
    select(
        CPM.cutoff, residual.model, AD.definition.column,
        z.score.counts, n.SVs, SVs.included,
        sample.strategy, 
        CRDID, detected.in.groups, PeakList
    ) %>%
    unnest(PeakList) %>%
    left_join(all.peak.differential.results.df, by=join_by(CPM.cutoff, residual.model, PeakID)) %>%
    nest(Peaks=c(PeakID, logFC, p.value, p.adjust))
all.differential.CRD.results.df <- 
    ALL_CRD_DIFFERENTIAL_STATS_TSV_FILEPATH %>%
    read_tsv(show_col_types=FALSE) %>% 
    filter(!is.na(p.value)) %>% 
    select(
        CPM.cutoff, residual.model,
        z.score.counts, n.SVs, SVs.included,
        CRD.sample.strategy, contrast, method,
        CRDID, p.value, p.adjust
    ) %>%
    mutate(contrast=str_split_i(contrast, '-', i=1)) %>% 
    dplyr::rename('sample.strategy'=CRD.sample.strategy)
CRDs.df <- 
    all.differential.CRD.results.df %>% 
    filter(method == 'deltaSLE') %>% 
    inner_join(
        all.merged.CRDs.df,
        by=intersect(colnames(all.differential.CRD.results.df), colnames(all.merged.CRDs.df))
    ) %>%
    mutate(
        omni.p.value=
            pmap_dbl(.l=., .f=function(Peaks, ...) { min(1, nrow(Peaks) * min(Peaks$p.value)) }),
        omni.p.adjust=
            pmap_dbl(.l=., .f=function(Peaks, ...) { min(1, nrow(Peaks) * min(Peaks$p.adjust)) })
    )

############################################################
# stageR
############################################################
# all.peak.differential.results.df
# all.merged.CRDs.df %>% dplyr::count(sample.strategy, meanClusterSize, detected.in.groups)
CRDs.df %>%
# all.differential.CRD.results.df
    nest(
        stager.input=
            -c(
                CPM.cutoff, residual.model, 
                z.score.counts, n.SVs, SVs.included,
                AD.definition.column, sample.strategy,
                contrast, contrast.levels, method,
            )
    ) %>% 
    {.} -> tmp; stager.input=tmp$stager.input[[1]]; stager.input
run_stageR <- function(
    stager.input,
    screen.metric='omni.p.value',
    alpha=0.05, ...){
    # n.sig.CRDs <- stager.input %>% filter(omni.p.value > alpha) %>% nrow()
    # stager.input %>% filter(omni.p.value < alpha) %>% nrow()
    # stager.input %>% filter(omni.p.adjust < alpha) %>% nrow()
    # stager.input %>% filter(p.value < alpha) %>% nrow()
    # stager.input %>% filter(p.adjust < alpha) %>% nrow()
    pScreen <- 
        stager.input %>% select(CRDID, !!sym(screen.metric)) %>%
        dplyr::rename('gene'=CRDID) %>% 
        as.data.frame() %>% column_to_rownames('gene')
    pConfirmation <- 
        stager.input %>% select(CRDID, Peaks) %>% unnest(Peaks) %>% 
        select(PeakID, p.value) %>%
        dplyr::rename('transcript'=PeakID) %>% 
        deframe()
    CRD.peak.mapping <- 
        stager.input %>% select(CRDID, Peaks) %>% unnest(Peaks) %>% 
        select(CRDID, PeakID) %>% distinct() %>% 
        dplyr::rename('gene'=CRDID, 'transcript'=PeakID) %>%
        as.data.frame()
    n.sig.CRDs <- 
        length(pScreen < alpha)
    if (n.sig.CRDs < 2) { return(NULL) }
    stager.results <- 
        stageR::stageRTx(
            pScreen = pScreen,
            pConfirmation = pConfirmation,
            pScreenAdjusted = FALSE,
            tx2gene=CRD.peak.mapping
        ) %>% 
        stageR::stageWiseAdjustment(
            object = .,
            method = "dte",
            alpha = alpha,
            allowNA = FALSE
        ) %>% 
        stageR::getAdjustedPValues(
            order = TRUE,
            onlySignificantGenes = FALSE
        )

    stager.results$CRD <- stager.results$geneID
    stager.results$PeakID <- stager.results$txID
}


