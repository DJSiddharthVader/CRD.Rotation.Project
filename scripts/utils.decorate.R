############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(GenomicRanges)
    library(plyranges)
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
        if (sample.strategy == 'All.Samples') {
            . %>% filter(!is.na(!!sym(AD.definition.column)))
        } else if (AD.definition.column == 'AD_CERAD_withDLB' & sample.strategy == 'Control.Samples') {
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
tidy_CRD_scores <- function(
    scores.obj,
    ...){
    scores.obj %>% 
    enframe(
        name='meanClusterSize',
        value='scores'
    ) %>% 
    unnest(scores) %>%
    as_tibble() %>%  
    dplyr::rename('ClusterID'='cluster', 'chr'='chrom') 
}

tidy_CRD_clusters <- function(
    clusters.obj,
    ...) {
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
    unnest(peak.clusters)
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
    separate_wider_delim(
        CRDID,
        delim='_',
        names=c('meanClusterSize', 'chr', NA)
    ) %>% 
    dplyr::rename('ClusterID'='cluster') %>% 
    arrange(meanClusterSize, seqnames, start, end)
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
    CRDID.prefix='',
    # CRDID.suffix='',
    ...) {
    # decorate.blob.filepath=tmp$decorate.blob.filepath[[1]]; CRDID.prefix=tmp$CRDID.prefix[[1]];
    results_dir <- dirname(decorate.blob.filepath)
    decorate.blob <- readRDS(decorate.blob.filepath)
    blobs.df <- 
        tribble(
            ~is.decorate.filtered, ~CRDs.blob,
            'all',                 decorate.blob$all.clusters,
            'filtered',            decorate.blob$filtered.clusters,
        )
    ##################################################################
    # Save mapping of individual PeakIDs <-> CRDIDs to a file 
    blobs.df %>% 
    mutate(
        CRDs=
            future_pmap(
                .l=.,
                CRDID.prefix=CRDID.prefix,
                .f=
                    function(CRDs.blob, CRDID.prefix, is.decorate.filtered, ...)
                    inner_join(
                        CRDs.blob %>% convert_PeakList_to_Locations(all.peak.locations.granges),
                        CRDs.blob %>% tidy_CRD_clusters() %>% nest(PeakList=c(PeakID)),
                        by=join_by(meanClusterSize, chr, ClusterID)
                    ) %>% 
                    {
                        if (is.decorate.filtered == 'all') {
                            inner_join(
                                .,
                                decorate.blob$all.LEFs %>% tidy_CRD_scores(),
                                by=join_by(meanClusterSize, chr, id, ClusterID)
                            )
                        } else {
                            .
                        }
                    } %>% 
                    mutate(CRDID=glue('{CRDID.prefix}{meanClusterSize}-{chr}-{start}-{end}'))
            )
    ) %>% 
    select(-c(CRDs.blob)) %>% 
    unnest(CRDs) %>% 
    select(-c(ClusterID, id)) %>% 
    nest(CRDs=-c(is.decorate.filtered, meanClusterSize)) %>% 
    mutate(
        results_dir=
            file.path(
                results_dir,
                glue('meanClusterSize_{meanClusterSize}'),
                glue('is.decorate.filtered_{is.decorate.filtered}')
            )
    ) %>% 
    future_pmap(
        .l=.,
        .f=
            function(results_dir, CRDs, is.decorate.filtered, ...){
                # CRDs=tmp2$CRDs[[1]]
                dir.create(results_dir, recursive=TRUE, showWarnings=FALSE)
                # Save list of constituent peaks for eahc 
                CRDs %>% 
                select(CRDID, PeakList) %>% 
                unnest(PeakList) %>%
                write_tsv(file.path(results_dir, 'PeakList-decorate.CRDs.tsv'))
                # Save genomic coordinates for each set of CRDs to a file
                CRDs %>% 
                select(CRDID, seqnames, start, end, width, strand, chr) %>% 
                write_tsv(file.path(results_dir, 'Locations-decorate.CRDs.tsv'))
                # Save CRD correlation summary stats
                if (is.decorate.filtered == 'all') {
                    CRDs %>% 
                    select(CRDID, N, LEF, mean_abs_corr, quantile75, quantile90, quantile95) %>%
                    write_tsv(file.path(results_dir, 'decorate.CRDs.LEFs.tsv'))
                } 
                invisible(NULL)
            }
    )
}

############################################################
# Merge CRDs and annotated with merging groups per merged CRD
############################################################
combine_CRD_Locations_and_PeakLists <- function(
    PeakList.filepath,
    Locations.filepath,
    ...){
    # PeakList.filepath=tmp$PeakList.filepath[[1]]; Locations.filepath=tmp$Locations.filepath[[1]];
    inner_join(
        PeakList.filepath %>% read_tsv(show_col_types=FALSE) %>% nest(PeakList=c(PeakID)),
        Locations.filepath %>% read_tsv(show_col_types=FALSE),
        by=join_by(chr, CRDID)
    )
}

merge_CRDs <- function(
    CRDs.to.merge,
    CRDID.prefix='',
    ...){
    # CRDs.to.merge=tmp$CRDs.to.merge[[1]]; CRDID.prefix=tmp$CRDID.prefix[[1]];
    merging_column <- setdiff(colnames(CRDs.to.merge), 'CRDs')[[1]]
    # For all CRDs across all values of the merging_column
    CRDs.to.merge %>% 
    unnest(CRDs) %>% 
    unnest(PeakList) %>%
    nest(PeakList=c(PeakID)) %>% 
    # Nest all other CRD features into a single column for handling by reduce_ranges
    select(-c(CRDID, detected.in.groups)) %>% 
    nest(
        CRD.data=
            -c(
                seqnames, start, end, width, strand,
                !!sym(merging_column)
            )
    ) %>% 
    as_granges() %>% 
    # merge all CRDs that overlap > 0 bp, concat all associated metadata from constituent CRDs
    # for each merged CRD
    reduce_ranges(
        detected.in.groups=unique(!!sym(merging_column)),
        CRD.data=c(CRD.data)
    ) %>% 
    as_tibble() %>% 
    mutate(
        # list of all unique CRD attributes of all CRDs merged to create this merged CRD (row)
        detected.in.groups=pmap_chr(.l=list(detected.in.groups), .f=~ paste0(sort(.x), collapse='+')),
        # Take all the peak lists from each CRD being merged and dedup the constituent peaks
        CRD.data=
            pmap(
                .l=list(CRD.data),
                .f=
                    ~ .x %>% 
                    bind_rows() %>% 
                    unnest(PeakList) %>% 
                    distinct() %>%
                    nest(PeakList=c(PeakID))
            ),
    ) %>% 
    unnest(CRD.data) %>% 
    mutate(ClusterID=row_number()) %>% 
    mutate(CRDID=glue('{CRDID.prefix}{chr}-{start}-{end}')) %>% 
    select(-c(ClusterID))
}

############################################################
# Compute features per CRD + summary stats over all Peaks per CRD
############################################################
compute_frequency <- function(vector){
    tibble(category.name=vector) %>%
    count(category.name) %>% 
    mutate(freq=n / sum(n)) %>% 
    pivot_longer(
        c(n, freq),
        names_to='metric',
        values_to='CRD.value'
    ) %>% 
    select(category.name, metric, CRD.value)
}

compute_CRD_stats <- function(
    Locations.filepath,
    PeakList.filepath,
    peak.metadata,
    differential.results,
    ...){
    # Locations.filepath=tmp$Locations.filepath[[1]]; PeakList.filepath=tmp$PeakList.filepath[[1]]; differential.results=tmp$differential.results[[1]]; peak.metadata=tmp$peak.metadata[[1]];
    CRDs <- 
        combine_CRD_Locations_and_PeakLists(
            PeakList.filepath,
            Locations.filepath
        )
    # compute summary stats over all peak features
    all.stats.df <- 
        CRDs %>% 
        unnest(PeakList) %>%
        left_join(
            peak.metadata,
            by=join_by(chr, PeakID)
        ) %>% 
        group_by(CRDID) %>%
        summarize(
            CRD.length=unique(width),
            CRD.nPeaks=n(),
            CRD.density=n() / unique(width),
            CRD.peak.coverage=sum(peak.length) / unique(width),
            freqs.peak.type=list(compute_frequency(peak.annotation.simple)),
            freqs.gene.type=list(compute_frequency(gene.type)),
            freqs.merge.groups=list(compute_frequency(detected.in.groups)),
            across(
                .cols=
                    c(
                        peak.length,
                        gene.length,
                        peak.distanceToTSS,
                    ),
                .names='{.col}-{.fn}',
                .fns=STAT_SUMMARY_FNCS_LIST
            )
        )
    # Per CRD, for all Peaks compute frequency of individual peak + linked-gene annotations
    peak.frequencies <- 
        all.stats.df %>% 
        select(CRDID, c(starts_with('freqs.'))) %>% 
        pivot_longer(
            starts_with('freqs.'),
            names_to='category.type',
            names_prefix='freqs.',
            values_to='freq.df'
        ) %>%
        unnest(freq.df) %>% 
        mutate(CRD.feature=glue('{category.type}-{category.name}-{metric}')) %>%
        select(CRDID, CRD.feature, CRD.value)
    # For each CRD + each contrast compute summary stats over 
    # differential statistics across all OCRs 
    differential.stats.df <- 
        CRDs %>% 
        unnest(PeakList) %>%
        left_join(
            differential.results,
            by=join_by(chr, PeakID)
        ) %>% 
        mutate(log10.p.adjust=-log10(p.adjust)) %>% 
        group_by(CRDID, contrast) %>%
        summarize(
            .groups='drop',
            across(
                .cols=c(logFC, log10.p.adjust),
                .names='{.col}-{.fn}',
                .fns=STAT_SUMMARY_FNCS_LIST
            )
        ) %>% 
        pivot_longer(-c(CRDID, contrast), names_to='metric', values_to='CRD.value') %>% 
        mutate(CRD.feature=glue('{contrast}-{metric}')) %>%
        select(CRDID, CRD.feature, CRD.value)
    # bind all summary stats together, including % of peaks of specific categories
    list(
        'CRD.stats'=
            all.stats.df %>% 
            select(CRDID, CRD.length, CRD.nPeaks, CRD.peak.coverage) %>%
            pivot_longer(-c(CRDID), names_to='CRD.feature', values_to='CRD.value'),
        'peak.summary.stats'=
            all.stats.df %>% 
            select(
                -c(
                    CRD.length,
                    CRD.nPeaks,
                    CRD.peak.coverage,
                    starts_with('freqs.')
                )
            ) %>%
            pivot_longer(-c(CRDID), names_to='CRD.feature', values_to='CRD.value'),
        'peak.frequencies'=peak.frequencies,
        'peak.differential.stats'=differential.stats.df 
    ) %>%
    enframe(
        name='CRD.feature.type',
        value='CRD.stats'
    ) %>%
    unnest(CRD.stats)
}

############################################################
# differetial CRD Analysis w/ decorate
############################################################
compute_differential_CRD_results <- function(
    peak.residuals.mx,
    CRDs,
    contrast.indicator,
    min.CRD.size,
    method.corr,
    test.methods,
    ...) {
    # For each CRD
    CRDs %>% 
    mutate(nPeaks=pmap_int(.l=list(PeakList), .f=nrow)) %>% 
    # Remvoe if CRD is too small
    filter(nPeaks >= min.CRD.size) %>% 
    # Try every test method for every CRD
    cross_join(tibble(method=test.methods)) %>% 
    # Compute CRD testing results i.e. compare peak abundance matrices between the two contrast levels
    # using the peaks specified by each CRD
    # turn results into tidy tibble
    mutate(
        test.results=
            future_pmap_dfr(
                .l=.,
                # .progress=TRUE,
                # only compare samples applicable to the specified contrast
                contrast.peak.residuals.mx=peak.residuals.mx[, rownames(contrast.indicator)],
                sample.groups=as.integer(contrast.indicator[[1]]),
                method.corr=method.corr,
                .f=
                    function(
                        PeakList,
                        method,
                        contrast.peak.residuals.mx,
                        sample.groups,
                        method.corr,
                        ...) {
                        corrMatrix.test( 
                            Y=t(contrast.peak.residuals.mx[PeakList$PeakID, ]),
                            method=method,
                            group=sample.groups,
                            method.corr=method.corr
                        ) %>%
                        enframe(name='differential.stat', value='differential.value') %>%
                        pivot_wider(names_from=differential.stat, values_from=differential.value) %>% 
                        unnest(where(is.list))
                    }
            )
    ) %>% 
    unnest(test.results) %>%
    dplyr::rename(
        # stat'=stat,
        # count'=count,
        # sign=sign,
        'p.value'=pVal
    ) %>%
    select(-c(PeakList))
}

compute_all_differential_CRD_results <- function(
    residuals.filepath, coords.filepath, all.sample.metadata,
    SVs.filepath, SVs.included,
    PeakList.filepath, Locations.filepath,
    contrast.variable, contrast.numerator, contrast.denominator,
    # test.methods=c('Delaneau', 'deltaSLE', 'Box'),
    ...){
    # paste0('row.index=1; ', paste0(colnames(tmp), '=tmp$', colnames(tmp), '[[row.index]]', collapse='; '))
    # Load ATAC data
    peak.data <- 
        subset_peak_data(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=all.sample.metadata,
            AD.definition.column='AD_CERAD_withDLB',
            sample.strategy='All.Samples'
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
    # Only keep samples with valid values for the contrast I cam comparing
    contrast.indicator <- 
        peak.data$sample.metadata %>%
        filter(!!sym(contrast.variable) %in% c(contrast.numerator, contrast.denominator)) %>%
        filter(SampleID %in% colnames(peak.residuals.mx)) %>% 
        mutate(
            !!contrast.variable := 
                factor(
                    !!sym(contrast.variable),
                    levels=c(contrast.denominator, contrast.numerator)
                )
        ) %>%
        select(all_of(c('SampleID', contrast.variable))) %>%
        as.data.frame() %>%
        column_to_rownames('SampleID')
    # Load CRDs to compare
    CRDs <- 
        combine_CRD_Locations_and_PeakLists(
            PeakList.filepath,
            Locations.filepath
        )
    # Now use custom wrapper to compute differential CRD stats with decorate
    compute_differential_CRD_results(
        peak.residuals.mx=peak.residuals.mx,
        CRDs=CRDs,
        contrast.indicator=contrast.indicator,
        ...
        # test.methods=test.methods,
        # min.CRD.size=min.CRD.size,
        # method.corr=method.corr,
    )
}

