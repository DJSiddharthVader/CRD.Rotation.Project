#############################
# Relevant directories
#############################
# PROJECT_DIR <- "/sc/arion/projects/Microglia/CRD.Analysis.Sid.Rotation"
RESIDUALS_DIR <- 
    here(
        '..',
        "Combined_RNAseq_ATACseq",
        "ATACseq_analyses",
        "Data_processing",
        "files"
    )
SAMPLE_METADATA_RDS_FILEPATH <- 
    file.path(RESIDUALS_DIR, "AllInfo_ATACseq_processing_209x362.RDs")
SAMPLE_METADATA_TSV_FILEPATH <- 
    here("all.sample.metadata.tsv")
DIFFERENTIAL_OCR_RESULTS_RDS_FILEPATH <- 
    here(
        '..',
        "Combined_RNAseq_ATACseq",
        "ATACseq_analyses",
        "Differential_analyses",
        "GCv38_ATACseq_over45y_3CPM",
        "DEG_results",
        "CombinedDEGresult_GCv38_ATACseq_over45y_3CPM_List_08_20_26.RDs"
    )
DIFFERENTIAL_OCR_RESULTS_TSV_FILEPATH <- 
    here(
        'results',
        'differential.OCR.results.tsv'
    )
#############################
# Results output dirs
#############################
SVA_RESULTS_DIR                  <- here("results", "peak.SVs")
ELBOW_RESULTS_DIR                <- here("results", "elbow.data")
SV_VARIANCEPARTITION_RESULTS_DIR <- here("results", "SV.variancePartitions")
SV_PEAK_STATS_DIR                <- here("results", "SV.peak.stats")
CRD_RESULTS_DIR                  <- here("results", "decorate.CRDs")
