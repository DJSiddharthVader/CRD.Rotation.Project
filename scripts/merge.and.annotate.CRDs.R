############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
all.sample.metadata <- 
    load_sample_metadata()
# differential OCR data
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
    nest(
        differential.results=
            -c(
                CPM.cutoff, sample.strategy, residual.model,
                contrast, contrast.levels,
            )
    ) %>%
    select(-c(sample.strategy))

############################################################
# Merge any CRDs that overlap, across meanClusterSizes
############################################################
# merge CRDs within each sample set across meanClustersizes
# so any ovelrapping CRDs (> 0 bp overlapping) we take the union of all of them
merged.CRDs.df <- 
    # List all CRD results files
    list_all_results_files_in_set(set.name='CRDs') %>% 
    # Tidily combine coordinates + peaklists for each CRD, 1 per row
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                .f=combine_CRD_Locations_and_PeakLists,
                .progress=TRUE
            )
    ) %>% 
    # combine all CRDs at different meanClusterSizes into a single tibble per input
    nest(CRDs.to.merge=c(meanClusterSize, CRDs)) %>% 
        head(1) %>% 
        # {.$CRDs.to.merge[[1]]} -> CRDs.to.merge
    # merge all overlapping CRDs across meanClusterSizes
    mutate(
        merged.CRDs=
            future_pmap(
                .l=.,
                .f=merge_overlapping_CRDs,
                merging_column='meanClusterSize',
                .progress=TRUE
            )
    ) %>%
    select(-c(CRDs.to.merge)) %>% 
    add_column(meanClusterSize='merged')

############################################################
# "Super Merge" any merged CRDs that overlap across sample subsets
############################################################
# So we take all the CRDs defined for the same input across all sample sets
# i.e. AD, Control, All samples take the union of any set of overlapping CRDs
# and keep all individual CRDs, so shared CRDs become the union across sample sets
# and all condition-specific CRDs are retained
super.merged.CRDs.df <- 
    merged.CRDs.df %>% 
    filter(sample.strategy != 'No.Others') %>% 
    # combine all CRDs from different input sample sets into a single tibble per input
    dplyr::rename('CRDs'=merged.CRDs) %>% 
    nest(CRDs.to.merge=c(sample.strategy, CRDs)) %>% 
    # This also annotated which sample sets the CRD was specifically detected in
    # i.e. All.Samples, Control.Samples, AD.Samples or any combination if the CRD
    # was detected in multiple sample sets
    # This is stored in the columns `detected.sample.strategy`
    mutate(
        merged.CRDs=
            future_pmap(
                .l=.,
                .f=merge_overlapping_CRDs,
                merging_column='sample.strategy',
                .progress=TRUE
            )
    ) %>%
    select(-c(CRDs.to.merge)) %>% 
    add_column(sample.strategy='merged')

############################################################
############################################################

############################################################
# Save all merged CRDs to appropriately nested files
############################################################
merged.CRDs.df
super.merged.CRDs.df
list(
    merged=merged.CRDs.df,
    super.merged=super.merged.CRDs.df
) %>%
    bind_rows() %>% 
    cross_join(tibble(format=c('PeakList', 'Locations'))) %>% 
    mutate(
        results_dir=
            file.path(
                CRD_RESULTS_DIR,
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
                glue('is.decorate.filtered_{is.decorate.filtered}')
            ),
        results_file=file.path(results_dir, glue('{format}-decorate.CRDs.tsv'))
    ) %>% 
        {.}
    future_pmap(
        .l=.,
        .f=
            function(results_file, ...) {
                merged.CRDs %>% 
                mutate(
                    ClusterID=nrow(),
                    CRDID=glue(CRDID_GLUE_DEFINITION)
                ) %>% 
                {
                    if (format == 'PeakList') {
                        select(.,
                            c(
                                meanClusterSize, 
                                chr,
                                ClusterID,
                                CRDID,
                                PeakID
                            )
                        )
                    } else if (format == 'Locations') {
                        select(.,
                            c(
                                seqnames
                                start
                                end
                                width
                                strand
                                meanClusterSize
                                chr
                                ClusterID
                                CRDID
                            )
                        )
                    }
                } %>% 
                write_tsv(results_file)
            },
        .progress=TRUE
)

############################################################
# Compute summary statistics over sets of CRDs
############################################################
plan(multisession, workers=TOTAL_CORES)
all.CRD.stats.df <- 
    check_cached_results(
        results_file=ALL_CRD_STATS_TSV_FILEPATH,
        force_redo=TRUE,
        results_fnc=
            function(){
                all.CRDs.df %>%
                inner_join(
                    all.peak.differential.results.df,
                    by=
                        join_by(
                            CPM.cutoff,
                            residual.model,
                        )
                ) %>% 
                mutate(
                    CRD.stats=
                        future_pmap(
                            .l=.,
                            .f=compute_CRD_stats,
                            .progress=TRUE
                        )
                ) %>%
                select(-c(differential.results, CRDs)) %>% 
                unnest(CRD.stats) %>% 
                filter(
                    CRD.feature %in% 
                        c(
                            'CRD.length',              
                            'CRD.nPeaks',              
                            # 'min-logFC',               
                            # 'q05-logFC',               
                            # 'q25-logFC',               
                            'mean-logFC',              
                            # 'med-logFC',               
                            # 'q75-logFC',               
                            # 'q95-logFC',               
                            # 'max-logFC',               
                            'var-logFC',               
                            # 'sd-logFC',                
                            'total-logFC',             
                            # 'min-log10.p.adjust',      
                            # 'q05-log10.p.adjust',      
                            # 'q25-log10.p.adjust',      
                            'mean-log10.p.adjust',     
                            # 'med-log10.p.adjust',      
                            # 'q75-log10.p.adjust',      
                            # 'q95-log10.p.adjust',      
                            # 'max-log10.p.adjust',      
                            # 'var-log10.p.adjust',      
                            # 'sd-log10.p.adjust',       
                            'total-log10.p.adjust',    
                            # 'min-peak.length',         
                            # 'q05-peak.length',         
                            # 'q25-peak.length',         
                            'mean-peak.length',        
                            # 'med-peak.length',         
                            # 'q75-peak.length',         
                            # 'q95-peak.length',         
                            # 'max-peak.length',         
                            # 'var-peak.length',         
                            # 'sd-peak.length',          
                            'total-peak.length',       
                            # 'min-peak.distanceToTSS',  
                            # 'q05-peak.distanceToTSS',  
                            # 'q25-peak.distanceToTSS',  
                            'mean-peak.distanceToTSS', 
                            # 'med-peak.distanceToTSS',  
                            # 'q75-peak.distanceToTSS',  
                            # 'q95-peak.distanceToTSS',  
                            # 'max-peak.distanceToTSS',  
                            # 'var-peak.distanceToTSS',  
                            # 'sd-peak.distanceToTSS',   
                            # 'total-peak.distanceToTSS',
                            # 'q05-gene.length',         
                            # 'min-gene.length',         
                            # 'q25-gene.length',         
                            'mean-gene.length'
                            # 'med-gene.length',         
                            # 'q75-gene.length',         
                            # 'q95-gene.length',         
                            # 'max-gene.length',         
                            # 'var-gene.length',         
                            # 'sd-gene.length',          
                            # 'total-gene.length',       
                        )
                )
        }
    )

