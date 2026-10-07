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
            df %>% mutate(!!col.name := factor(!!sym(col.name), levels=col.levels))
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

list_all_ATAC_datasets <- function(){
    # First list residual ATAC-seq count matrices that can be used to annotate CRDs
    list_all_ATAC_residual_sets() %>% 
    # get data specifying peak coordinates for each set of residuals
    left_join(
        list_all_ATAC_coordinates_sets(),
        by=join_by(CPM.cutoff)
    )
}

############################################################
# Listing all results files + metadta prefedined named sets of results 
############################################################
list_all_results_files_in_set <- function(
    set.name,
    ...){
    if (set.name == 'SVs') {
        SVA_RESULTS_DIR %>% 
        parse_results_filelist(
            ...,
            suffix='peak.residual.SVs.tsv'
        ) %>%
        dplyr::rename('SVs.filepath'=filepath) %>% 
        select(-c(filename))
    } else if (set.name == 'SVs.elbow') {
        ELBOW_RESULTS_DIR %>% 
        # CRD_RESULTS_DIR %>% 
        parse_results_filelist(
            ...,
            suffix='elbow.cluster.data.tsv',
            filename.column.name='filename'
        ) # %>% select(-c(filename))
    } else if (set.name == 'SVs.variancePartition') {
        SV_VARIANCEPARTITION_RESULTS_DIR %>% 
        parse_results_filelist(
            ...,
            suffix='SV.variance.partition.results.tsv',
            filename.column.name='filename'
        ) %>% 
        select(-c(filename))
    } else if (set.name == 'CRD.blobs') {
        CRD_RESULTS_DIR %>%
        parse_results_filelist(
            ...,
            suffix='decorate.cluster.blob.rds'
        ) %>% 
        select(-c(filename))
    # } else if (set.name == '') {
    } else {
        stop(glue('Invalid results set name: {set.name}'))
    }
}

############################################################
# Load differential abundance results for individual OCRs
############################################################
load_differential_OCR_results <- function(...){
    check_cached_results(
        ...,
        results_file=DIFFERENTIAL_OCR_RESULTS_TSV_FILEPATH,
        results_fnc=
            function() {
                DIFFERENTIAL_OCR_RESULTS_RDS_FILEPATH %>%
                readRDS() %>%
                    {.} -> tmp; tmp
                enframe(
                    name='model.name',
                    value='model'
                ) %>%
                mutate(contrast=pmap(.l=list(model), .f=~ .x$ALL %>% names() %>% list())) %>%
                unnest(contrast) %>% 
                mutate(contrast.type=ifelse(str_detect(contrast, '_for_'), 'discrete', 'continuous')) %>% 
                separate_wider_delim(
                    contrast,
                    delim='_for_',
                    names=c(NA, 'contrast.variable'),
                    too_few='align_start',
                    cols_remove=FALSE
                ) %>% 
                mutate(
                    differential.results.df=
                        pmap(
                            .l=.,
                            .f=
                                function(model, variable.names, ...){
                                    model$ALL[[variable.names]] %>% 
                                    as_tibble() %>%
                                    dplyr::rename(
                                        'chr'=chromosome_name,
                                        'EnsemblID'=ensembl_gene_id,
										"peak.start"=start_position,
										"peak.end"=end_position,
                                        "peak.strand"=strand,
										"peak.annotation"=annotation,
                                        "peak.annotation.simple"=annotationSimple,
                                        "gene.chr"=geneChr,
										"gene.start"=geneStart,
										"gene.end"=geneEnd,
                                        "gene.length"=geneLength,
										"gene.strand"=geneStrand,
										"gene.name"=Associated.Gene.Name,
										"p.value"=P.Value,
                                        "p.adjust"=adj.P.Val

                                    ) %>% 
                                    mutate(peak.length=peak.end-peak.start)
                                }
                        )
                )
            }
    ) %>% 
    filter(contrast.type == 'discrete') %>% 
    add_column(sample.strategy='All.Samples') %>% 
    select(
        PeakID,
        EnsemblID, gene.name,
        chr, 
        peak.start,
        peak.end,
        # peak.strand,
        peak.annotation,
        peak.annotation.simple,
        # gene.chr,
        gene.start,
        gene.end,
        # gene.length,
        # gene.strand,
        # AveExpr, z.std, t,
        logFC,
        p.value, p.adjust
    )
}

