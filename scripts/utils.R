############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(tictoc)
    library(glue)
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
    show_col_types=FALSE,
    silence=FALSE,
    results_fnc,
    ...){
    # Now check if the results file exists and load it
    tic()
    if (is.null(results_file)) {
        if (!silence) { message("No results file, will just return data") }
        results <- results_fnc(...)
        return_data <- TRUE
    } else {
        output.filetype <- results_file %>% str_extract('\\.[^\\.]*$') 
        if (!silence) { message(output.filetype) }
        if (output.filetype == '.rds') {
            load_fnc <- readRDS
            save_fnc <- saveRDS
        } else if (output.filetype %in% c('.txt', '.tsv')) {
            load_fnc <- partial(read_tsv, show_col_types=show_col_types)
            save_fnc <- write_tsv
        } else {
            stop(glue('Invalid file extesion: {output.filetype}'))
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
            results <- results_fnc(...) %T>% save_fnc(results_file)
            # If results dont exist or force_redo is TRUE compute + cache results
            # Assumes save_fnc is of fomr save_fnc(result_object, filename)
        }
    }
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
            separate_longer_delim(
                .,
                fileinfo,
                delim='/'
            ) %>%
            mutate(
                fileinfo=
                    ifelse(
                        grepl(suffix_pattern, fileinfo),
                        paste(
                            filename.column.name,
                            fileinfo,
                            sep=param_delim
                        ) %>% 
                        str_remove(suffix),
                        fileinfo
                    )
            ) %>% 
            separate_wider_delim(
                fileinfo,
                delim='_',
                too_many="merge",  # in case filenames have delim inside
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
    ) %>% 
    mutate(
        age=as.integer(age),
        CDR_3levels=
            factor(
                CDR_3levels,
                levels=c('Healthy', 'MCI', 'Dementia')
            ),
        Braak_3levels=
            factor(
                Braak_3levels,
                levels=c('Braak_Max2', 'Braak_3to4', 'Braak_5plus')
            ),
        AD_CERAD_withDLB=
            factor(
                AD_CERAD_withDLB,
                levels=c('Control', 'AD', 'Other')
            )
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

