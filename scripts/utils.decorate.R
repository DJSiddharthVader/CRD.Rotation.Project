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
        load_peak_locations()
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

clean_decorate_param_columns <- function(df){
    decorate.params.list <- colnames(DECORATE_HYPER_PARAMS_DF) %>% {.[. != 'meanClusterSize']}
    {
        if (length(setdiff(decorate.params.list, colnames(df))) == 0) {
            df %>% 
            mutate(decorate.params=glue(paste0('{', decorate.params.list, '}', collapse='+'))) %>% 
            select(-all_of(decorate.params.list))
        } else {
            df
        }
    } %>% 
    mutate(meanClusterSize=factor(meanClusterSize))
}

############################################################
# Generate decorate CRDs
############################################################
tidy_CRD_scores <- function(scores.obj){
    scores.obj %>% 
    enframe(
        name='meanClusterSize',
        value='scores'
    ) %>% 
    unnest(scores) %>%
    as_tibble() # %>% select(-c(id))
}

tidy_CRD_clusters <- function(clusters.obj) {
    # every row is assignment for a single peak to a CRD
    # columns: meanClusterSize, chr, PeakID, CRDID
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
    mutate(CRDID=glue(CRDID_GLUE_DEFINITION))
    # mutate(CRDID=glue('{chr}.{CRDID}'))
}

convert_PeakList_to_Locations <- function(
    treeListClusters,
    all.peak.locations.granges,
    ...){
    # use decorate function to do this automatically
    getClusterRanges(all.peak.locations.granges, treeListClusters) %>% 
    # convert to a tibble for saving the file 
    names_to_column('CRDID') %>% 
    as_tibble() %>%
    select(-c(id, cluster)) %>% 
    separate_wider_delim(
        CRDID,
        delim='_',
        names=c('meanClusterSize', 'chr', 'ClusterID')
    ) %>% 
    mutate(CRDID=glue(CRDID_GLUE_DEFINITION)) %>% 
    arrange(meanClusterSize, seqnames, start, end, CRDID)
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
    # Get "scores" for each CRD from analyzing the OCR correlation matrix 
    # defines the internal correlation and "strength" of each CRD
    tree.scores <- 
        tree.list %>% 
        scoreClusters(
            treeList=.,
            treeListClusters=all.tree.list.clusters,
            BPPARAM=BPPARAM
        )
    # Remove clusters scores below thresholds and collapse any clusters that overlap 
    filtered.tree.list.clusters <- 
        tree.scores %>% 
        retainClusters( 
            metric=filterMetric,
            cutoff=filterMetricCutoff
        ) %>% 
        # get retained clusters
        filterClusters(all.tree.list.clusters, .) %>% 
        # collapse redundant clusters i.e. any pair of clusters with > jaccardCutoff
        collapseClusters(
            peak.locations,
            jaccardCutoff=jaccardCutoff
        )
    
    # Score the filtered + colllapsed CRDs 
    # THINK THIS BREAKS BECAUSE OF A BUG IN DECORATE SEE ./scratch.R
    # filtered.scores <- 
    #     tree.list %>% 
    #     scoreClusters(
    #         treeList=.,
    #         treeListClusters=filtered.tree.list.clusters,
    #         BPPARAM=BPPARAM
    #     )
    # return all data to be saved as an .rds object
    # some of these can saved as a individual tsvs, which we do later
    list(
        all.clusters=all.tree.list.clusters,
        all.LEFs=tree.scores,
        filtered.clusters=filtered.tree.list.clusters,
        # filtered.LEFs=filtered.scores,
        tree.list=tree.list
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
            ...
            # # colnames(tmp) %>% tibble(name=.)
            # residuals.filepath=residuals.filepath,
            # residual.model  =residual.model    ,
            # CPM.cutoff      =CPM.cutoff        ,
            # coords.filepath =coords.filepath   ,
            # sample.strategy =sample.strategy   ,
            # AD.definition.colum=AD.definition.colum,
            # SVs.filepath    =SVs.filepath      ,
            # n.SVs           =n.SVs             ,
            # z.score.counts  =z.score.counts    ,
            # SVs.included    =SVs.included      ,
            # adjacentCount   =adjacentCount     ,
            # method.corr     =method.corr       ,
            # clusterMethod   =clusterMethod     ,
            # meanClusterSize =meanClusterSize   ,
            # filterMetric    =filterMetric      ,
            # filterMetricCutoff=filterMetricCutoff,
            # jaccardCutoff   =jaccardCutoff     
        )
    if (!is.null(p)) { p() }
    message(names(results.obj))
    return(results.obj)
}

unpack_decorate_blob <- function(
    decorate.blob.filepath,
    all.peak.locations.granges,
    ...) {
    # decorate.blob.filepath=tmp$decorate.blob.filepath[[1]];
    results_dir <- dirname(decorate.blob.filepath)
    decorate.blob <- readRDS(decorate.blob.filepath)
    # peak <-> CRD (cluster) mapping
    decorate.blob$all.clusters %>%
        tidy_CRD_clusters() %>%
        write_tsv(file.path(results_dir, 'all-PeakList-decorate.CRDs.tsv'))
    # CRDs decinfed by their genomic coordinates
    decorate.blob$all.clusters %>%
        convert_PeakList_to_Locations(all.peak.locations.granges) %>% 
        write_tsv(file.path(results_dir, 'all-Locations-decorate.CRDs.tsv'))
    # CRD self-correlation stats 
    decorate.blob$all.LEFs %>%
        tidy_CRD_scores() %>%
        write_tsv(file.path(results_dir, 'all-decorate.CRD.LEFs.tsv'))
    # post filtering+collapsing by decorate
    decorate.blob$filtered.clusters %>%
        tidy_CRD_clusters() %>%
        write_tsv(file.path(results_dir, 'filtered-PeakList-decorate.CRDs.tsv'))
    # CRDs decinfed by their genomic coordinates
    decorate.blob$filtered.clusters %>%
        convert_PeakList_to_Locations(all.peak.locations.granges) %>% 
        write_tsv(file.path(results_dir, 'filtered-Locations-decorate.CRDs.tsv'))
    # stats for filtered CRDs
    decorate.blob$filtered.LEFs %>%
        tidy_CRD_scores() %>%
        write_tsv(file.path(results_dir, 'filtered-decorate.CRD.LEFs.tsv'))
    invisible(NULL)
    # decorate.blob$tree.list %>%
    #     tidy_decorate_tree.list() %>% 
    #     write_tsv(file.path(results_dir, 'decorate.tree.list.tsv'))
}

############################################################
# Post-process CRDs (merge + annotate)
############################################################
combine_CRD_Locations_and_PeakLists <- function(
    PeakList.filepath,
    Locations.filepath,
    ...){
    inner_join(
        PeakList.filepath %>% 
            read_tsv(show_col_types=FALSE) %>% 
            nest(PeakList=c(PeakID)),
        Locations.filepath %>% 
            read_tsv(show_col_types=FALSE),
        by=
            join_by(
                meanClusterSize, 
                chr,
                ClusterID,
                CRDID
            )
    ) %>%
    nest(CRDs=-c(meanClusterSize))
}

load_all_CRD_data <- function(
    filter_fnc=NULL,
    ...){
    # List all sets of CRDs estimated
    check_cached_results(
        results_file=ALL_CRD_PEAKLISTS_AND_LOCATIONS_FILEPATH,
        ...,
        filter_fnc=filter_fnc,
        results_fnc=
            function(filter_fnc, ...) {
                list_all_results_files_in_set(set.name='CRDs') %>%
                {
                    if (!is.null(filter_fnc)) {
                        filter_fnc(.)
                    } else {
                        .
                    }
                } %>% 
                # Load PeakLists and coordinates separately and 
                # then map together for downstream analysis (1 CRD per row)
                mutate(
                    CRDs=
                        future_pmap(
                            .l=.,
                            .f=combine_CRD_Locations_and_PeakLists,
                            .progress=TRUE
                        )
                ) %>% 
                unnest(CRDs)
            }
    ) %>% 
    clean_decorate_param_columns() %>% 
    # clean up columns
    select(
        CPM.cutoff, residual.model, AD.definition.column, sample.strategy,
        z.score.counts, n.SVs, SVs.included,
        decorate.params, meanClusterSize, is.filtered,
        CRDs
    )
}

merge_overlapping_CRDs <- function(
    CRDs,
    all.peak.metadata.df,
    ...){
    # row.index=1; CRDs=all.CRDs.df$CRDs[[row.index]];
    CRDs %>% 
    as_granges() %>% 
    # reduce_ranges()
    reduce_ranges(
        chr=unique(chr),
        ClusterIDs=list(c(ClusterID)),
        CRDIDs=list(c(CRDID)),
        PeakList=list(c(PeakList))
    ) %>% 
        {.}
}

compute_frequency <- function(vector, category.name){
    tibble(category=vector) %>%
    count(category) %>% 
    mutate(category=glue('{category.name}-{category}')) %>% 
    muatate(freq=n / n()) %>% 
    list()
}

compute_CRD_stats <- function(
    CRDs,
    differential.results,
    ...){
    CRDs %>% 
    unnest(PeakList) %>%
    left_join(
        differential.results,
        # copy=TRUE,
        by=
            join_by(
                chr,
                PeakID
            )
    ) %>% 
    mutate(log10.p.adjust=-log10(p.adjust)) %>% 
    group_by(CRDID) %>%
    summarize(
        CRD.length=unique(width),
        CRD.nPeaks=n(),
        # peak.type.freqs=compute_frequency(peak.type, 'peak.type'),
        # gene.type.freqs=compute_frequency(gene.type, 'gene.type'),
        across(
            .cols=
                c(
                    logFC,
                    log10.p.adjust,
                    peak.length,
                    peak.distanceToTSS,
                    gene.length
                ),
            .names='{.fn}-{.col}',
            .fns=STAT_SUMMARY_FNCS_LIST
        )
    ) %>% 
    unnest(ends_with('.freqs')) %>% 
    pivot_longer(
        -c(CRDID),
        names_to='CRD.feature',
        values_to='CRD.value',
    )
}

