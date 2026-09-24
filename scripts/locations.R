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
SAMPLE_METADATA_RDS_FILEPATH <- file.path(RESIDUALS_DIR, "AllInfo_ATACseq_processing_209x362.RDs")
SAMPLE_METADATA_TSV_FILEPATH <- here("all.sample.metadata.tsv")
DECORATE_RESULTS_DIR         <- here("results", "decorate")
SVA_RESULTS_DIR              <- here("results", "peak.SVs")
ELBOW_RESULTS_DIR            <- here("results", "elbow.data")
CRD_RESULTS_DIR              <- here("results", "decorate.CRDs")
