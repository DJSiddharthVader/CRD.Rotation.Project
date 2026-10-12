############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(tictoc)
    library(glue)
    library(forcats)
    # library(optparse)
    # library(furrr)
    # library(plyranges)
})

############################################################
# Utilities
############################################################
check_cached_results <- function(
    results_file,
    force_redo=FALSE,
    return_data=TRUE,
    silence=FALSE,
    p=NULL,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        if (silence) { 
            results <- invisible(results_fnc(...))
        } else {
            message("No results file, will just return data") 
            results <- results_fnc(...)
        }
        return_data <- TRUE
    } else {
        output.filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (!silence) { message(output.filetype) }
        if (output.filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output.filetype %in% c('.txt', '.tsv')) {
            load_fnc <- read_tsv
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {output.filetype}'))
        }
        # silence functions is specified
        if (silence) {
            load_fnc <- quietly(load_fnc)
            save_fnc <- quietly(save_fnc)
        }
        if (file.exists(results_file) & !force_redo) {
            if (!silence) { message(glue('{results_file} exists, not recomputing results')) }
            if (return_data) {
                if (!silence) { message('Loading cached results...') }
                results <- load_fnc(results_file)
            }
        # or force recompute+cache the data and return it
        } else {
            if (file.exists(results_file) & force_redo) {
                if (!silence) { message(glue('{results_file} exists, recomputing results anyways')) }
            } else {
                if (!silence) { message(glue('No cached results, generating: {results_file}')) }
            }
            dir.create(dirname(results_file), recursive=TRUE, showWarnings=FALSE)
            if (silence) { 
                # results <- invisible(results_fnc(...))
                results <- invisible(results_fnc(...)) %T>% save_fnc(results_file)
            } else {
                results <- results_fnc(...) %T>% save_fnc(results_file)
            }
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
    if (!is.null(p)) { p() } # update progress bar if specified
    # Or dont save the results, just return the results
    if (!silence) { toc() }
    if (return_data) {
        return(results)
    } else { 
        return(invisible(NULL))
    }
}

parse_results_filelist <- function(
    input_dir,
    suffix,
    pattern=NA,
    param_delim='_',
    filename.column.name='filename',
    parse_filepath_to_columns=TRUE,
    ...){
    # !!NOTICE!!
    # This will break if any parameter_dir name has a param_delim character in the name or value, 
    # This shouldnt  break if the filename has a single param_delim character in it 
    # i.e. min_resolution-0.45 is fine when param_delim='-'
    #      min_resolution_0.45 is will break with any value for param_delim not '_'
    suffix_pattern <- glue('*{suffix}$')
    # List all results files that exist
    input_dir %>% 
    list.files(
        pattern=suffix_pattern,
        recursive=TRUE,
        full.names=FALSE
    ) %>% 
    tibble(fileinfo=.) %>%
    mutate(filepath=file.path(input_dir, fileinfo)) %>% 
    {
        if (is.na(pattern)) {
            .
        } else {
            filter(., grepl(pattern, filepath))
        }
    } %>% 
    # Extract param info from directory names into structured columns
    {
        if (parse_filepath_to_columns) {
            mutate(., !!filename.column.name := basename(str_remove(fileinfo, suffix))) %>% 
            mutate(fileinfo=dirname(fileinfo)) %>% 
            separate_longer_delim(
                .,
                fileinfo,
                delim='/'
            ) %>%
            separate_wider_delim(
                fileinfo,
                delim=param_delim,
                too_many="merge",  # in case filenames have delim inside
                too_few="align_start",
                names=
                    c(
                        'Parameter',
                        'Value'
                    )
            ) %>%
            pivot_wider(
                names_from=Parameter,
                values_from=Value
            )
        } else {
            .
        }
    } %>% 
    # Fix column types based on content
    readr::type_convert()
}

set_col_levels <- function(df, col.name) {
    if (col.name  %in% colnames(df)) {
        col.levels <- ALL_VARIABLE_ORDERINGS[[col.name]]
        if (!is.null(col.levels)) {
            df %>% mutate(!!col.name := fct_drop(factor(!!sym(col.name), levels=col.levels)))
        } else {
            df %>% mutate(!!col.name := as.integer(!!sym(col.name)))
        }
    } else {
        df
    }
}

order_variable_for_plotting <- function(
    df,
    to.skip=NULL,
    to.keep=NULL,
    ...) {
    columns.to.transform <- 
        colnames(df) %>%
        intersect(names(ALL_VARIABLE_ORDERINGS)) %>% 
        {
            if (!is.null(to.skip)) {
                setdiff(., to.skip)
            } else {
                .
            }
        } %>% 
        {
            if (!is.null(to.keep)) {
                intersect(., to.skip)
            } else {
                .
            }
        }
    # print(columns.to.transform)
    for (col.name in columns.to.transform) {
        df <- df %>% set_col_levels(col.name=col.name)
    } 
    return(df)
}

############################################################
# Loading ATAC Data
############################################################
load_sample_metadata <- function(...){
    check_cached_results(
        ...,
        results_file=SAMPLE_METADATA_TSV_FILEPATH,
        results_fnc=
            function(){
                SAMPLE_METADATA_RDS_FILEPATH %>% 
                readRDS() %>% 
                rownames_to_column('SampleID') %>% 
                as_tibble() %>% 
                select(
                    SampleID,
                    # Sample_ID,
                    Person_ID,
                    sex,
                    age,
                    RACE,
                    race_color,
                    Biobank,
                    readCountInfo_Total,
                    finalReadCount,
                    lib.size,
                    ctcf_fos,
                    mergingDesigns,
                    PrimCauseDeath,
                    PMI_mins,
                    Sepsis_Dx,
                    Secondary_Dx,
                    Parkinson_Dx,
                    # DLB_Dx,
                    # Dementia_MCI_CDR,
                    # CDR_Combined,
                    # Braak_Score,
                    # CERAD_Combined,
                    ApoE_score,
                    ApoE4count,
                    AD_CERAD_withDLB,
                    Braak_3levels,
                    CDR_3levels,
                    Final.Dx
                )
            }
    )
}

list_all_ATAC_residual_sets <- function(){
    RESIDUALS_DIR %>%
    list.files(full.names=TRUE) %>%
    tibble(residuals.filepath=.) %>%
    mutate(info=residuals.filepath %>% basename()) %>% 
    filter(str_detect(info, '^Resid_.*.RDs')) %>% 
    mutate(info=info %>% str_remove('.RDs$')) %>% 
    separate_wider_delim(
        info,
        delim='_',
        names=c(NA, 'residual.model', NA, NA, 'CPM.cutoff', NA, NA),
        too_few='align_start',
        too_many='merge'
    ) %>%
    filter(str_detect(CPM.cutoff, 'CPM'))
}

list_all_ATAC_coordinates_sets <- function(){
    RESIDUALS_DIR %>%
    list.files(full.names=TRUE) %>%
    tibble(coords.filepath=.) %>%
    mutate(info=coords.filepath %>% basename()) %>% 
    filter(str_detect(info, '^expObjNonLow_ATACseq_.*.RDs')) %>% 
    mutate(info=info %>% str_remove('.RDs$')) %>% 
    separate_wider_delim(
        info,
        delim='_',
        names=c(NA, NA, 'dims', 'CPM.cutoff'),
        too_few='align_start',
        too_many='merge'
    ) %>%
    filter(str_detect(CPM.cutoff, 'CPM')) %>%
    select(-c(dims))
}

load_peak_locations <- function(coords.filepath, ...) {
    coords.filepath %>% 
    readRDS() %>% 
    {.$genes} %>% 
    as.data.frame() %>% 
    makeGRangesFromDataFrame()
}

list_all_ATAC_datasets <- function(){
    # First list residual ATAC-seq count matrices that can be used to annotate CRDs
    list_all_ATAC_residual_sets() %>% 
    # get data specifying peak coordinates for each set of residuals
    left_join(
        list_all_ATAC_coordinates_sets(),
        by=join_by(CPM.cutoff)
    )
}

load_all_peak_metadata <- function(
    CPM.cutoff='3CPM',
    ...){
    check_cached_results(
        ...,
        results_file=PEAK_METADATA_TSV_FILEPATH,
        results_fnc=
            function(CPM.cutoff) {
                list_all_ATAC_coordinates_sets() %>% 
                mutate(
                    peaks=
                        future_pmap(
                            .l=list(coords.filepath),
                            .f= ~ .x %>% readRDS() %>% {.$genes} %>% as_tibble(),
                            .progres=TRUE
                        )
                ) %>%
                unnest(peaks) %>% 
                mutate(peak.length=end-start) %>% 
                dplyr::rename(
                    'chr'=seqnames,
                    "peak.start"=start,
                    "peak.end"=end,
                    "peak.strand"=strand,
                    "peak.distanceToTSS"=distanceToTSS,
                    "peak.annotation"=annotation,
                    "peak.annotation.simple"=annotationSimple,
                    'gene.EnsemblID'=geneId,
                    'transcript.EnsemblID'=transcriptId,
                    "gene.name"=Associated.Gene.Name,
                    # "gene.chr"=geneChr,
                    "gene.start"=geneStart,
                    "gene.end"=geneEnd,
                    "gene.length"=geneLength,
                    "gene.strand"=geneStrand,
                    'gene.type'=Gene.type,
                    'description'=Description
                )
            }
    ) %>% 
    filter(CPM.cutoff == CPM.cutoff) %>% 
    select(
        PeakID,
        gene.EnsemblID, gene.name, 
        # transcript.EnsemblID,
        chr, peak.start, peak.end, peak.length,
        # peak.strand,
        peak.annotation, peak.annotation.simple,
        peak.distanceToTSS,
        # gene.chr,
        gene.start, gene.end, gene.length,
        # gene.strand,
        gene.type
    )
}

load_all_differential_OCR_results <- function(...){
    check_cached_results(
        ...,
        results_file=DIFFERENTIAL_OCR_TSV_FILEPATH,
        results_fnc=
            function() {
                DIFFERENTIAL_OCR_RDS_FILEPATH %>%
                readRDS() %>%
                enframe(
                    name='differential.model',
                    value='model'
                ) %>%
                mutate(
                    residual.model=
                        case_when(
                            differential.model == 'AgeModel'   ~ 'KeepDxAge',
                            differential.model == 'NoAgeModel' ~ 'KeepDx',
                        )
                ) %>% 
                unnest(model) %>% 
                mutate(
                    model=
                        future_pmap(
                            .l=list(model),
                            .f=~ .x %>% enframe(name='contrast', value='diff.results')
                        )
                ) %>% 
                unnest(model) %>% 
                separate_wider_delim(
                    contrast,
                    delim='_for_',
                    names=c('contrast.levels', 'contrast'),
                    # cols_remove=FALSE
                    too_few='align_end'
                ) %>% 
                mutate(contrast.type=ifelse(is.na(contrast.levels), 'continuous', 'discrete')) %>% 
                unnest(diff.results) %>% 
                dplyr::rename(
                    'chr'=chromosome_name,
                    "peak.start"=start_position,
                    "peak.end"=end_position,
                    "peak.strand"=strand,
                    "peak.annotation"=annotation,
                    "peak.annotation.simple"=annotationSimple,
                    "peak.distanceToTSS"=distanceToTSS,
                    "transcript.EnsemblID"=transcriptId,
                    "gene.EnsemblID"=ensembl_gene_id,
                    "gene.chr"=geneChr,
                    "gene.start"=geneStart,
                    "gene.end"=geneEnd,
                    "gene.length"=geneLength,
                    "gene.strand"=geneStrand,
                    "gene.name"=Associated.Gene.Name,
                    # "gene.type"=Gene.type,
                    "description"=Description,
                    "p.value"=P.Value,
                    "p.adjust"=adj.P.Val
                ) %>% 
                mutate(peak.length=peak.end-peak.start)
            }
    ) %>% 
    filter(contrast.type == 'discrete') %>% 
    add_column(
        CPM.cutoff='3CPM',
        sample.strategy='All.Samples'
    ) %>% 
    select(
        CPM.cutoff, sample.strategy, residual.model,
        PeakID,
        gene.EnsemblID, gene.name, # gene.type,
        chr, peak.start, peak.end, peak.length,
        # peak.strand,
        peak.annotation,
        peak.annotation.simple,
        peak.distanceToTSS,
        # description,
        # gene.chr,
        gene.start, gene.end, gene.length,
        # gene.strand,
        # AveExpr, z.std, t, B,
        # contrast.type,
        contrast, contrast.levels,
        logFC,
        p.value, p.adjust
    )
}

############################################################
# Listing all results files + metadta prefedined named sets of results 
############################################################
list_all_results_files_in_set <- function(
    set.name,
    ...){
    # List all filepaths to SVs estimated
    if (set.name == 'SVs') {
        SVA_RESULTS_DIR %>% 
        parse_results_filelist(suffix='peak.residual.SVs.tsv') %>%
        dplyr::rename('SVs.filepath'=filepath) %>% 
        select(-c(filename)) 
    # List computed variance of each peak's abundance across SVs residualized
    } else if (set.name == 'peak.variance.stats') {
        SV_PEAK_STATS_DIR %>% 
        parse_results_filelist(suffix='peak.stats.with.SVs.removed.tsv') %>% 
        dplyr::rename('peak.stats.filepath'=filepath) %>% 
        select(-c(filename))
    # List all filepaths to variancePartition results after iteratively regressing out SVs
    } else if (set.name == 'SVs.variancePartition') {
        SV_VARIANCEPARTITION_RESULTS_DIR %>% 
        parse_results_filelist(suffix='SV.variance.partition.results.tsv') %>% 
        dplyr::rename('SV.variancePartition.filepath'=filepath) %>%  
        select(-c(filename))
    # List all LEFs results extracted from blobs after iteratively regressing out SVs
    } else if (set.name == 'SVs.elbow') {
        CRD_RESULTS_DIR %>% 
        parse_results_filelist(suffix='decorate.CRD.LEFs.tsv') %>% 
        dplyr::rename('elbow.filepath'=filepath) %>% 
        select(-c(filename))
    # List all CRD objects with clusters + filtered clusters + LEFs
    } else if (set.name == 'CRD.blobs') {
        CRD_RESULTS_DIR %>%
        parse_results_filelist(suffix='decorate.blob.rds') %>% 
        dplyr::rename('decorate.blob.filepath'=filepath)  %>% 
        select(-c(filename))
    # List all CRDs results files with separate files for PeakLists and genomic positions
    } else if (set.name == 'CRDs') {
        CRD_RESULTS_DIR %>% 
        parse_results_filelist(suffix='-decorate.CRDs.tsv', filename.column.name='info') %>% 
        separate_wider_delim(info, delim='-', names=c('is.filtered', 'format')) %>% 
        mutate(format=glue('{format}.filepath')) %>% 
        pivot_wider(names_from=format, values_from=filepath)
    # } else if (set.name == '') {
    } else {
        stop(glue('Invalid results set name: {set.name}'))
    }
}

############################################################
# Load all CRDs, specified as CRDID <-> PeakID mappings
############################################################
load_all_CRDs_as_Peak_lists <- function() {
    decorate.params.list <- DECORATE_HYPER_PARAMS_DF %>% colnames() %>% {.[. != 'meanClusterSize']}
    list_all_results_files_in_set(set.name='CRD.PeakLists') %>% 
    # only load CRDs with chosen SVs regressed out 
    # filter(CPM.cutoff == '3CPM') %>% 
    # filter(z.score.counts) %>% 
    # filter(n.SVs == NUM_SVS_TO_GENERATE) %>% 
    # filter(SVs.included == NUM_SVS_TO_REMOVE) %>% 
    # filter(sample.strategy != 'No.Others') %>% 
    mutate(decorate.params=glue(paste0('{', decorate.params.list, '}', collapse='+'))) %>% 
    mutate(
        CRDs=
            future_pmap(
                .l=list(CRDs.filepath),
                .f=~ read_tsv(.x, show_col_types=FALSE),
                .progress=TRUE
            )
    ) %>%
    unnest(CRDs) %>%
    nest(CRDs=c(chr, PeakID, CRDID)) %>% 
    arrange(
        CPM.cutoff,
        residual.model,
        decorate.params, meanClusterSize,
        AD.definition.column, sample.strategy,
        is.decorate.filtered
    ) %>% 
    select(
        CPM.cutoff, residual.model, AD.definition.column, sample.strategy,
        z.score.counts, n.SVs, SVs.included,
        decorate.params, meanClusterSize, is.decorate.filtered,
        CRDs
    )
}

