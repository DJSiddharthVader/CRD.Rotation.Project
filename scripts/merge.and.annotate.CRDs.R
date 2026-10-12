############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'basic.imports.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
all.sample.metadata <- 
    load_sample_metadata()
# only keep results subset for testing
example_filter <- function(df){
    df %>% 
    filter(sample.strategy != 'No.Others') %>% 
    filter(CPM.cutoff == '3CPM') %>% 
    filter(z.score.counts) %>% 
    filter(n.SVs == NUM_SVS_TO_GENERATE) %>% 
    filter(SVs.included == NUM_SVS_TO_REMOVE)
    # select(-all_of(decorate.params.list)) %>% 
    # select(-c(CPM.cutoff, z.score.counts, n.SVs, SVs.included))
}

############################################################
# Load all generated CRD + Peak data
############################################################
plan(multisession, workers=TOTAL_CORES)
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
# Load CRDs across all conditions, including list of peaks and coordiantes for each CRD (1 per row)
all.CRDs.df <- 
    load_all_CRD_data(filter_fnc=example_filter)
    # load_all_CRD_data(force_redo=TRUE, filter_fnc=example_filter)

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

############################################################
