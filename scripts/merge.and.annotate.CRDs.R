############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))

############################################################
# Load all individually called CRDs post LEF filtering
############################################################
# Load All CRDs, 23s w/ 20 cores
plan(multisession, workers=TOTAL_CORES)
all.CRDs.df <- 
    # List all CRD results files
    list_all_results_files_in_set(set.name='CRDs') %>% 
    filter(sample.strategy %in% c('AD.Samples', 'Control.Samples', 'All.Samples')) %>% 
    filter(SVs.included == NUM_SVS_TO_REMOVE) %>% 
    # merging only makes sense for filterd clusters
    # otherwise all clusters willl just be the max meanClusterSize CRDs
    filter(is.decorate.filtered == 'filtered') %>% 
    # Tidily combine coordinates + peaklists for each CRD, 1 per row
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                .progress=TRUE,
                .f=
                    function(PeakList.filepath, Locations.filepath, ...) {
                        combine_CRD_Locations_and_PeakLists(
                            PeakList.filepath,
                            Locations.filepath
                        ) %>%
                        add_column(detected.in.groups=NA)
                    }
            )
    ) %>% 
    select(-c(ends_with('.filepath')))

############################################################
# Merge any CRDs that overlap, across meanClusterSizes
############################################################
# merge CRDs within each sample set across meanClustersizes
# so any ovelrapping CRDs (> 0 bp overlapping) we take the union of all of them
plan(multisession, workers=TOTAL_CORES)
merged.CRDs.df <- 
    all.CRDs.df %>% 
    nest(CRDs.to.merge=c(meanClusterSize, CRDs)) %>% 
    add_column(meanClusterSize='merged') %>% 
    mutate(CRDID.prefix=glue('{AD.definition.column}-{sample.strategy}-{meanClusterSize}-')) %>% 
    # merge all overlapping CRDs across meanClusterSizes
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                .f=merge_CRDs,
                .progress=TRUE
            )
    ) %>%
    select(-c(CRDs.to.merge, CRDID.prefix))

############################################################
# "Super Merge" any merged CRDs that overlap across sample subsets
############################################################
# So we take all the CRDs defined for the same input across all sample sets
# i.e. AD, Control, All samples take the union of any set of overlapping CRDs
# and keep all individual CRDs, so shared CRDs become the union across sample sets
# and all condition-specific CRDs are retained
plan(multisession, workers=TOTAL_CORES)
super.merged.CRDs.df <- 
    merged.CRDs.df %>%
    # combine all CRDs from different input sample sets into a single tibble per input
    nest(CRDs.to.merge=c(sample.strategy, CRDs)) %>% 
    add_column(sample.strategy='merged') %>% 
    mutate(CRDID.prefix=glue('{AD.definition.column}-{sample.strategy}-{meanClusterSize}-')) %>% 
    # This also annotates which sample sets the final "super-merged" CRD was specifically 
    # detected in # i.e. All.Samples, Control.Samples, AD.Samples or any combination 
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                .f=merge_CRDs,
                .progress=TRUE
            )
    ) %>%
    select(-c(CRDs.to.merge, CRDID.prefix))
# 
############################################################
# Save all merged CRDs to appropriately nested files
############################################################
tribble(
    ~ merging.level, ~CRDs,
    'separate',      all.CRDs.df %>% mutate(meanClusterSize=as.character(meanClusterSize)),
    'size.merged',   merged.CRDs.df,
    'sample.merged', super.merged.CRDs.df,
    ) %>%
    unnest(CRDs) %>% 
    cross_join(tibble(format=c('PeakList', 'Locations'))) %>% 
    mutate(
        results_dir=
            file.path(
                MERGED_CRD_RESULTS_DIR,
                glue('CPM.cutoff_{CPM.cutoff}'),
                glue('residual.model_{residual.model}'),
                glue('AD.definition.column_{AD.definition.column}'),
                glue('sample.strategy_{sample.strategy}'),
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
                glue('merging.level_{merging.level}')
            ),
        results_file=file.path(results_dir, glue('{format}-decorate.CRDs.tsv'))
    ) %>% 
    pivot_wider(
        names_from=format,
        names_prefix='results_file.',
        values_from=results_file
    ) %>% 
    future_pmap(
        .l=.,
        .progress=TRUE,
        .f=
            function(results_file.PeakList, results_file.Locations, results_dir, CRDs, ...) {
                # save CRDs as mapping of Peaks <-> CRDs 
                # for merged CRDs also keep the sub.CRDID they came from when merging
                dir.create(results_dir, showWarnings=FALSE, recursive=TRUE)
                CRDs %>% 
                select(
                    c(
                        chr,
                        CRDID,
                        PeakList
                    )
                ) %>%
                unnest(PeakList) %>%
                write_tsv(results_file.PeakList)
                # Save CRD as genomic coodrinates
                CRDs %>% 
                select(
                    c(
                        seqnames, start, end, width, strand,
                        CRDID,
                        chr,
                        detected.in.groups
                    )
                ) %>%
                write_tsv(results_file.Locations)
                invisible(NULL)
            }
    )

############################################################
# Compute summary statistics + features for each CRD across all sets/conditions
############################################################
# Load OCR metadata
all.peak.metadata.df <- 
    load_all_peak_metadata() %>%
    nest(peak.metadata=-c(CPM.cutoff))
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
        CPM.cutoff, sample.strategy, residual.model,
        PeakID, 
        chr, peak.start, peak.end, peak.length, 
        contrast, contrast.levels,
        logFC, p.value, p.adjust
    ) %>% 
    nest(differential.results=-c(CPM.cutoff, sample.strategy, residual.model))
# Compute summary stats over OCRs per CRD
plan(multisession, workers=TOTAL_CORES)
list_all_results_files_in_set(set.name='merged.CRDs') %>% 
    inner_join(
        all.peak.metadata.df,
        by=join_by(CPM.cutoff)
    ) %>% 
    inner_join(
        all.peak.differential.results.df,
        by=join_by(CPM.cutoff, residual.model)
    ) %>% 
    mutate(results_file=file.path(dirname(PeakList.filepath), 'CRD.stats.tsv')) %>% 
    future_pmap(
        .l=.,
        .f=check_cached_results,
        # force_redo=TRUE,
        return_data=FALSE,
        results_fnc=compute_CRD_stats,
        .progress=TRUE
    )

