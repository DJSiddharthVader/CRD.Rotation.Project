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

#############################
# Data Loading
#############################
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
        # Set read/write functions based on filetype
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
                    Sample_ID,
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
    tibble(filepath=.) %>%
    mutate(info=filepath %>% basename()) %>% 
    filter(str_detect(info, '^Resid_.*.RDs')) %>% 
    separate_wider_delim(
        info,
        delim='_',
        names=c(NA, 'model', NA, NA, 'CPM.cutoff', NA, NA),
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

############################################################
}

