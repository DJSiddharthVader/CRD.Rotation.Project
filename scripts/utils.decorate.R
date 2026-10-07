############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(decorate)
    library(limma)
    source(here('scripts', 'utils.SVs.R'))
})

############################################################
# Input PreProcessing
############################################################
select_sample_subset <- function(
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    ...){
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='AD.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='Control.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='No.Others'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='All.Samples'
    # AD.definition.column='AD_CERAD_withDLB'; sample.strategy='Make.Error'
    # decide which specific samples to select  based on inclusin criteria
    filter_cmd <- 
        if        (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'Control.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Control')
        } else if (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'AD.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'AD')
        } else if (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'No.Others') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) != 'Other')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'Control.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Braak_3to4')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'AD.Samples') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) == 'Braak5plus')
        } else if (AD.definition.column == 'Bradk_3levels' & sample.strategy == 'No.Others') {
            .  %>% dplyr::filter(!!sym(AD.definition.column) != 'BraakMax2')
        } else if (sample.strategy == 'All.Samples') {
            . %>% filter(!is.na(!!sym(AD.definition.column)))
        } else {
            stop(glue('sample inclusion criteria not implemented. you used {sample.strategy} to filter on the column {AD.definition.column}'))
        }
    # Now pull the Sample_IDs from the metadata for the selected samples
    sample.metadata %>%
    filter_cmd() %>% 
    pull(SampleID)
}

subset_peak_data <- function(
    residuals.filepath,
    coords.filepath,
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    ...){
    # load coordinates of all ATAC features
    peak.locations.gr <- 
        coords.filepath %>% 
        readRDS() %>% 
        {.$genes} %>% 
        as.data.frame() %>% 
        makeGRangesFromDataFrame()
        # select(seqnames, start, end)
    # load residual ATAC matrix
    peak.residuals.mx <- 
        residuals.filepath %>% 
        readRDS()
    # Only keep common peaks
    peaks.to.keep <- 
        intersect(
            names(peak.locations.gr),
            rownames(peak.residuals.mx)
        )
    # select which samples to include when annotating CRDs 
    samples.to.keep <- 
        select_sample_subset(
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # Now make sure rows are compatible
    # length(samples.to.keep); length(peaks.to.keep)
    # peak.locations.gr <- peak.locations.gr[peaks.to.keep, ]
    # peak.residuals.mx <- peak.residuals.mx[peaks.to.keep, samples.to.keep]
    list(
        sample.metadata=sample.metadata %>% filter(SampleID %in% samples.to.keep),
        residuals=peak.residuals.mx[peaks.to.keep, samples.to.keep],
        locations=peak.locations.gr[peaks.to.keep, ]
    )
}

############################################################
# Generate decorate clusters
############################################################
tidy_CRD_scores <- function(scores.obj){
    scores.obj %>% 
    enframe(
        name='meanClusterSize',
        value='scores'
    ) %>% 
    unnest(scores) %>%
    as_tibble() %>%
    select(-c(id))
}

tidy_CRD_clusters <- function(clusters.obj) {
    # every row is assignment for a single peak to a CRD
    # columns: meanClusterSize, chr, PeakID, ClusterID
    clusters.obj %>%
    names() %>% 
    tibble(meanClusterSize=.) %>% 
    mutate(
        peak.clusters=
            future_pmap(
                .l=.,
                .f=
                    function(meanClusterSize, ...) {
                        all.genome.clusters <- clusters.obj[[meanClusterSize]]
                        all.genome.clusters %>% 
                        names() %>% 
                        tibble(chr=.) %>% 
                        mutate(
                            clusters=
                                pmap(
                                    .l=.,
                                    .f=
                                        function(chr, ...){
                                            all.genome.clusters[[chr]] %>% 
                                            enframe(
                                                name='PeakID',
                                                value='ClusterID'
                                            )
                                        }
                                )
                        ) %>%
                        unnest(clusters)
                    }
            )
    ) %>%
    unnest(peak.clusters) %>% 
    mutate(ClusterID=glue('{chr}.{ClusterID}'))
}

generate_peak_cluster_with_decorate <- function(
    peak.residuals.mx,
    peak.locations,
    adjacentCount,
    method.corr,
    clusterMethod,
    meanClusterSize,
    jaccardCutoff,
    filterMetric,
    filterMetricCutoff,
    BPPARAM=SerialParam(),
    ...){
    # peak.residuals.mx=peak.residuals.mx; peak.locations=peak.data$locations; BPPARAM=SnowParam(N_CORES_PER_PROCESS)
    # Evaluate hierarchical clustering
    # adjacentCount is the number of adjacent peaks considered in correlation
    # use Spearman correlation to reduce the effects of outliers
    # peak.residuals.mx=peak.data$residuals; peak.locations=peak.data$locations;
    tree.list <- 
        peak.residuals.mx %>% 
        runOrderedClusteringGenome( 
            peak.locations,
            adjacentCount=adjacentCount,
            method.corr=method.corr
        )
    # Choose cutoffs and return clusters using multiple values for meanClusterSize 
    # Clusters corresponding to each parameter value are returned and then processed downstream
    # By using multiple parameter values, epigenetic features are included in clusters 
    # at different resolutions
    all.tree.list.clusters <- 
        tree.list %>% 
        createClusters(
            method=clusterMethod,
            meanClusterSize=meanClusterSize
        )
    tree.scores <- 
        tree.list %>% 
        scoreClusters(
            all.tree.list.clusters,
            BPPARAM=BPPARAM
        )
    filtered.tree.list.clusters <- 
        tree.scores %>% 
        retainClusters( 
            metric=filterMetric,
            cutoff=filterMetricCutoff
        ) %>% 
        # get retained clusters
        filterClusters(all.tree.list.clusters, .) %>% 
        # collapse redundant clusters
        collapseClusters(
            peak.locations,
            jaccardCutoff=jaccardCutoff
        )
    filtered.scores <- 
        tree.list %>% 
        scoreClusters(
            filtered.tree.list.clusters,
            BPPARAM=BPPARAM
        )
    # return all data
    list(
        tree.list=tree.list,
        all.clusters=all.tree.list.clusters,
        all.LEFs=tree.scores,
        filtered.clusters=filtered.tree.list.clusters,
        filtered.LEFs=filtered.scores
    )
}

run_decorate_pipeline <- function(
    residuals.filepath,
    coords.filepath,
    SVs.filepath,
    all.sample.metadata,
    AD.definition.column='AD_CERAD_withDLB',
    sample.strategy='All',
    SVs.included=Inf,
    # LEFs.only=FALSE,
    p=NULL,
    ...) {
    # paste0(c('row.index=1', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; ')), collapse='; ')
    # row.index=1; residuals.filepath=tmp$residuals.filepath[[row.index]]; residual.model=tmp$residual.model[[row.index]]; CPM.cutoff=tmp$CPM.cutoff[[row.index]]; coords.filepath=tmp$coords.filepath[[row.index]]; sample.strategy=tmp$sample.strategy[[row.index]]; AD.definition.column=tmp$AD.definition.column[[row.index]]; SVs.filepath=tmp$SVs.filepath[[row.index]]; n.SVs=tmp$n.SVs[[row.index]]; z.score.counts=tmp$z.score.counts[[row.index]]; SVs.included=tmp$SVs.included[[row.index]]; adjacentCount=tmp$adjacentCount[[row.index]]; method.corr=tmp$method.corr[[row.index]]; clusterMethod=tmp$clusterMethod[[row.index]]; meanClusterSize=tmp$meanClusterSize[[row.index]]; filterMetric=tmp$filterMetric[[row.index]]; filterMetricCutoff=tmp$filterMetricCutoff[[row.index]]; jaccardCutoff=tmp$jaccardCutoff[[row.index]]; results_dir=tmp$results_dir[[row.index]]; results_file=tmp$results_file[[row.index]]
    # subset + order data as defined
    peak.data <- 
        subset_peak_data(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=all.sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # Load SVs to residualize out 
    if (SVs.included > 0) {
        SVs.df <- 
            SVs.filepath %>%
            read_tsv(show_col_types=FALSE) %>%
            as.data.frame() %>%
            column_to_rownames('SampleID')
        peak.residuals.mx <-
            residualize_SVs(
                SVs.df=SVs.df,
                residuals.mx=peak.data$residuals,
                SVs.included=SVs.included
            )
    } else  if (SVs.included == 0) {
        peak.residuals.mx <- 
            peak.data$residuals
    } else {
        stop(glue('SVs.included must be an integer, passed: {SVs.included}'))
    }
    results.obj <- 
        generate_peak_cluster_with_decorate(
            peak.residuals.mx=peak.residuals.mx,
            peak.locations=peak.data$locations,
            # ...
            # colnames(tmp) %>% tibble(name=.)
            residuals.filepath=residuals.filepath,
            residual.model  =residual.model    ,
            CPM.cutoff      =CPM.cutoff        ,
            coords.filepath =coords.filepath   ,
            sample.strategy =sample.strategy   ,
            AD.definition.colum=AD.definition.colum,
            SVs.filepath    =SVs.filepath      ,
            n.SVs           =n.SVs             ,
            z.score.counts  =z.score.counts    ,
            SVs.included    =SVs.included      ,
            adjacentCount   =adjacentCount     ,
            method.corr     =method.corr       ,
            clusterMethod   =clusterMethod     ,
            meanClusterSize =meanClusterSize   ,
            filterMetric    =filterMetric      ,
            filterMetricCutoff=filterMetricCutoff,
            jaccardCutoff   =jaccardCutoff     
        )
    if (!is.null(p)) { p() }
    message(names(results.obj))
    return(results.obj)
}

############################################################
# Misc posterity code
############################################################
sva_docs_example_correlations <- function() {
    library(bladderbatch)
    data(bladderdata)
    pheno = pData(bladderEset)
    edata = exprs(bladderEset)
    mod = model.matrix(~as.factor(cancer), data=pheno)
    mod0 = model.matrix(~1,data=pheno)
    n.sv = num.sv(edata,mod,method="be"); n.sv
    svobj = sva(edata,mod,mod0,n.sv=n.sv)
    colnames(svobj$sv) <- paste0('SV', 1:n.sv)
    canCorPairs(
        paste(c('~ cancer', colnames(svobj$sv)), collapse='+'),
        data=bind_cols(svobj$sv, pheno)
    )
    n.sv = num.sv(edata,mod,method="leek"); n.sv
    svobj = sva(edata,mod,mod0,n.sv=n.sv)
    colnames(svobj$sv) <- paste0('SV', 1:n.sv)
    canCorPairs(
        paste(c('~ cancer', colnames(svobj$sv)), collapse='+'),
        data=bind_cols(svobj$sv, pheno)
    )
    # my function shich should be the same
    run_sva(
        sample.metadata=as_tibble(pheno ),
        counts.matrix=edata,
        full.model.vars=c('cancer'),
        reduced.model.vars=NULL,
    ) %>% 
    as.data.frame() %>% 
    column_to_rownames('SampleID') %>% 
    {
        canCorPairs(
            paste(c('~ cancer', colnames(.)), collapse='+'),
            data=bind_cols(., pheno)
        )
    }
}

