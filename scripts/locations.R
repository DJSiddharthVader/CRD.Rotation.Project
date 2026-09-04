#############################
# Relevant directories
#############################
# PROJECT_DIR <- "/sc/arion/projects/Microglia/CRD.Analysis.Sid.Rotation"
MICROGLIA_DIR <- "/sc/arion/projects/Microglia"
RESIDUALS_DIR <- 
    file.path(
        MICROGLIA_DIR,
        "Combined_RNAseq_ATACseq",
        "ATACseq_analyses",
        "Data_processing", 
        "files"
    )
SAMPLE_METADATA_RDS_FILEPATH <- 
    file.path(
        RESIDUALS_DIR,
        "AllInfo_ATACseq_processing_209x362.RDs"
    )

SAMPLE_METADATA_TSV_FILEPATH <- 
    file.path(
        PROJECT_DIR,
        "all.sample.metadata.tsv"
    )

