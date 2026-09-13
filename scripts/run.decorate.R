############################################################
# Dependencies
############################################################
library(here)
source(here('scripts', 'locations.R'))
source(here('scripts', 'utils.R'))
source(here('scripts', 'utils.decorate.R'))
# load patient metadata
sample.metadata <- load_sample_metadata()


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
