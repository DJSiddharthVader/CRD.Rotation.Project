############################################################
# Constants
############################################################
NUM_SVS_TO_GENERATE <- 14 # number of OCRs / number of samples
NUM_SVS_TO_REMOVE <- 5 # decided manually using elbow plots of SVs
AD_PHENOTYPE_COLUMNS <- 
    c(
        'AD_CERAD_withDLB', 
        'Braak_3levels',
        'CDR_3levels'
    )
ALL_PHENOTYPE_COLUMNS <- 
    c(
        AD_PHENOTYPE_COLUMNS,
        'age'
    )
RELEVANT_METADATA_COLUMNS <- 
    c(
        AD_PHENOTYPE_COLUMNS,
        'age',
        'PMI_mins',
        'Biobank',
        'sex',
        'RACE'
    )
# how to order these variable as factor for plotting (left to right in plots)
# NULL means convert to integer
ALL_VARIABLE_ORDERINGS <- 
    list(
        age=NULL,
        CDR_3levels=
            c(
                'Healthy',
                'MCI',
                'Dementia'
            ),
        Braak_3levels=
            c(
                'Braak_Max2',
                'Braak_3to4',
                'Braak_5plus'
            ),
        AD_CERAD_withDLB=
            c(
                'AD',
                'Control',
                'Other'
            ),
        residual.model=
            c(
                'KeepDxAge',
                'KeepDx',
                'Full'
            ),
        sample.strategy=
            c(
                'merged',
                'All.Samples',
                'AD.Samples',
                'Control.Samples',
                'No.Others'
            ),
        meanClusterSize=
            c(
                'merged',
                '100',
                '80',
                '50',
                '25',
                '10'
            ),
        detected.in.groups=
            c(
                "10",
                "10+25",                                 
                "10+25+50",                             
                "10+25+50+80",                          
                "10+25+50+80+100",                      
                "10+50",                                
                "10+80",                                
                "10+100",                               
                "10+25+80",                             
                "10+50+80",                             
                "10+25+100",                            
                "10+50+100",                            
                "10+80+100",                            
                "10+25+50+100",                         
                "10+25+80+100",                         
                "AD.Samples",                           
                "Control.Samples",                      
                "All.Samples",                          
                "AD.Samples+All.Samples+Control.Samples",
                "AD.Samples+All.Samples",               
                "All.Samples+Control.Samples",          
                "AD.Samples+Control.Samples"
            )
    )

############################################################
# Hyper-param combinations
############################################################
RESIDUAL_MATRIX_SAMPLE_PARAMS_DF <- 
    tibble(
        # Different inputs to generate CRD annotations from
        sample.strategy=
            c(
                'Control.Samples',   # only annotated CRDs using Control samples
                'AD.Samples',        # only annotate CRDs using AD samples
                'No.Others',         # All Control + AD samples
                'All.Samples'        # All Control + AD + Other samples
            ),
        # which metadata column we use to define sample conditions
        AD.definition.column=
            c(
                # 'CDR_3levels',     # Dementia, Control, MCI
                # 'Braak_3levels',   # Braak_5plus, Braak_3to4, Braak_Max2
                'AD_CERAD_withDLB' # AD, Control, Other
            )
    )
DECORATE_HYPER_PARAMS_DF <- 
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
DIFFERENTIAL_DECORATE_HYPER_PARAMS_DF <- 
    tribble(
        ~contrast.variable, ~contrast.denominator,  ~contrast.numerator,
        'AD_CERAD_withDLB', 'Control',              'AD', 
        'Braak_3levels',    'Braak_Max2',           'Braak_3to4',
        'Braak_3levels',    'Braak_3to4',           'Braak_5plus',
        'Braak_3levels',    'Braak_Max2',           'Braak_5plus',
        'CDR_3levels',      'Healthy',              'MCI',
        'CDR_3levels',      'Healthy',              'Dementia',
        'CDR_3levels',      'MCI',                  'Dementia',
        'Sex',              'Female',               'Male',
        'Biobank',          'MSBB',                 'Rush',
    ) %>% 
    mutate(contrast=glue('{contrast.variable}-{contrast.numerator}-{contrast.denominator}'))

############################################################
# Common functions/summary stats
############################################################
STAT_SUMMARY_FNCS_LIST <- 
    c(
        # 'min'=min, 
        # 'q05'=~ quantile(.x, na.rm=TRUE, probs=0.05),
        # 'q25'=~ quantile(.x, na.rm=TRUE, probs=0.25),
        'mean'=mean,
        # 'med'=median, 
        # 'q75'=~ quantile(.x, na.rm=TRUE, probs=0.75),
        # 'q95'=~ quantile(.x, na.rm=TRUE, probs=0.95),
        # 'max'=max,
        # 'var'=var,
        'sd'=sd,
        'total'=~ sum(abs(.x), na.rm=TRUE)
    )

