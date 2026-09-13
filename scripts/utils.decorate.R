############################################################
# Dependencies
############################################################
suppressPackageStartupMessages({
    library(sva)
    library(GenomicRanges)
    library(decorate)
    library(variancePartition)
})

############################################################
# Run SVA on peak residuals
############################################################
run_sva <- function(
    sample.metadata,
    counts.matrix,
    full.model.vars,
    reduced.model.vars=NULL,
    ...) {
    # Run SVA to identify SVs and return SV matrix
    # First convert variable lists to formula objs
    reduced.fml <- 
        if (is.null(reduced.model.vars)) {
            formula('~ 1')
        } else {
            reduced.model.vars %>%
            paste(collapse='+') %>%
            sprintf('~%s', .)
            formula() 
        }
    # effects  + covaraites
    full.fml <- 
        full.model.vars %>%
        paste(collapse='+') %>%
        sprintf('~%s', .) %>% 
        formula() 
    # Run SVA 
    # full.fml %>% model.matrix(sample.metadata) %>% dim()
    # reduced.fml %>% model.matrix(sample.metadata) %>% dim()
    svs <- 
        tryCatch(
            {
                counts.matrix %>%
                as.matrix() %>%
                {.[rowSums(.) > 0, ]} %>%
                # svaseq(
                sva(
                    # full model matrix
                    full.fml %>% model.matrix(sample.metadata),
                    # covariate only model matrix
                    reduced.fml %>% model.matrix(sample.metadata),
                    # ...
                ) %>%
                {.$sv} %>%
                set_colnames(paste0('SV', 1:ncol(.))) %>% 
                as.data.frame() %>% 
                as_tibble()
            },
            error=function(cond){
                message('Estimating SVs failed with following error:')
                message(conditionMessage(cond))
            }
        )
}

run_sva_on_peak_residuals <- function(
    residuals.filepath,
    coords.filepath,
    sample.metadata,
    AD.definition.column,
    sample.strategy,
    full.SV.model.vars=c('AD_CERAD_withDLB', 'Braak_3levels', 'CDR_3levels', 'age'),
    reduced.SV.model.vars=NULL,
    ...) {
    # subset + order data as defined
    peak.data <- 
        pick_samples_to_use(
            residuals.filepath=residuals.filepath,
            coords.filepath=coords.filepath,
            sample.metadata=sample.metadata,
            AD.definition.column=AD.definition.column,
            sample.strategy=sample.strategy
        )
    # remove any rows with NAs in any model variable
    sample.metadata <- 
        sample.metadata %>%
        filter(SampleID %in% colnames(peak.data$residuals)) %>%
        filter(
            if_all(
                all_of(full.SV.model.vars),
                ~ !is.na(.)
            )
        )
    # remove any variables that are uniform across all samples
    for (model.var in full.SV.model.vars) {
        if (length(unique(sample.metadata[[model.var]])) == 1) {
            # sample.metadata <- sample.metadata %>% select(-c(!!sym(model.var)))
            full.SV.model.vars <- full.SV.model.vars[full.SV.model.vars != model.var]
            message(glue('removed {model.var} from SV model since it is singular'))
        }
    }
    # estimate SVs from peak residuals 
    run_sva(
        # counts.matrix=peak.data$residuals[1:1000, sample.metadata$SampleID],
        counts.matrix=peak.data$residuals[, sample.metadata$SampleID],
        sample.metadata=sample.metadata,
        full.model=full.SV.model.vars,
        reduced.model=reduced.SV.model.vars
    )
}

list_all_SVs <- function(){
    SVA_RESULTS_DIR %>% 
    parse_results_filelist(
        suffix='peak.residual.SVs.tsv',
        filename.column.name='filename'
    ) %>%
    select(-c(filename))
}

adjust_residuals_with_svs <- function(
