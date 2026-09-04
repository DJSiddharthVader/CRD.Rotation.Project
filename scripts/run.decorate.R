#############################
# Dependencies
#############################
PROJECT_DIR <- here()
source(file.path(PROJECT_DIR, 'scripts', 'locations.R'))

#############################
# Run decorate
#############################
all.residual.sets.df <- 
    list_all_ATAC_residual_sets()

all.residual.sets.df %>%
    pmap(
        .l=.,
        .f=run_decorate,
        resids.matrix,
        resids.granges,
        adjacentCount=500,
        method.corr="spearman",
        clusterMethod="meanClusterSize",
        meanClusterSize=c(20, 30, 40, 50),
        filterMetric="LEF",
        filterMetricCutoff=0.15,
        jaccardCutoff=0.9
    ) 
