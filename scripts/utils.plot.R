###################################################
# Dependencies
###################################################
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(ggh4x)
    library(ggridges)
    library(GGally)
    library(scales)
    library(ggragged)
    # library(ggExtra)
    library(ggside)
    # library(ggpointdensity)
    # library(viridis)
    library(ComplexUpset)
    library(quarto)
})

###################################################
# Handling/formatting plots
###################################################
make_ggtheme <- function(...){
    theme(
        panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.background=element_blank(),
        axis.line=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks=
            element_line(
                color="black",
                linewidth=1/2.13
            ),
        axis.ticks.length=unit(3, "pt"),
        axis.title=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=10
            ),
        axis.text=
            element_text(
                family="sans",
                face="bold",
                color="black",
                size=8
            ),
        ...
    )
}

build_axis_fnc <- function(
    figure,
    scale.mode,
    scale.axis,
    limits){
    # fetch axis data to check if data is discerete or not
    scale.data <- rlang::eval_tidy(rlang::quo_squash(figure@mapping[[scale.axis]]), figure@data)
    # Determine axis type
    axis.type <- 
        case_when(
            scale.mode == 'pct'      ~ 'continuous', # percentage
            scale.mode == 'm'        ~ 'continuous', # large numbers
            scale.mode == 'mb'       ~ 'continuous', # larger numbers in megabases
            scale.mode == 'log10'    ~ 'log10',
            is.character(scale.data) ~ 'discrete',
            is.factor(scale.data)    ~ 'discrete',
            scale.mode == 'discrete' ~ 'discrete',
            # scale.mode == 'binned'   ~ 'binned',
            scale.mode == ''         ~ 'continuous'
        )
    # build axis scaling functiona name for ggplot2 and eval to return function object
    glue('scale_{scale.axis}_{axis.type}') %>%
    as.symbol() %>%
    eval()

}

scale_axis <- function(
    figure,
    scale.axis='x',
    scale.mode='',
    log.base=10,
    axis.label.accuracy=0.1,
    n.breaks=NULL,
    limits=NULL,
    expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    # figure out which axis function to call
    axis_fnc <- 
        build_axis_fnc(
            figure=figure,
            scale.axis=scale.axis,
            scale.mode=scale.mode,
            limits=limits
        )
    # Set axis scaling/labeling based on scale.mode argument
    # Scaling for percentages
    if (scale.mode == 'pct') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=label_percent(),
            ...
        )
    # Human scaling for large numbers e.g. 1K, 1M, 1B
    } else if (scale.mode == 'm') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=
                label_number(
                    scale_cut=cut_short_scale(),
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    # Human scaling in bytes e.g. 1Kb, 1Mb, 1Gb. 
    # Useful for genomic coordinates + sizes 
    } else if (scale.mode == 'mb') {
        figure +
        axis_fnc(
            expand=expand,
            n.breaks=n.breaks,
            limits=limits,
            labels=
                label_bytes(
                    units="auto_si",
                    accuracy=axis.label.accuracy
                ),
            ...
        )
    # log scale + include minor ticks, uses base 10 by default
    } else if (scale.mode == 'log') {
        figure +
        axis_fnc(
            expand=expand,
            limits=limits,
            guide='axis_logticks',
            labels=
                label_log(
                    base=log.base,
                    digits=max(1, -log10(axis.label.accuracy)),
                    signed=FALSE
                ),
            ...
        )
    # Discrete set of labeld i.e. DEL,WT
    } else if (scale.mode == 'discrete') {
        figure +
        axis_fnc(expand=expand)
    # No scaling, just set limits and label rounding
    } else if (scale.mode == '') {
        if (is.null(limits)) {
            figure + 
            axis_fnc(
                expand=expand,
                labels=
                    function(x) {
                        format(x, digits=max(1, -log10(axis.label.accuracy)))
                    }
            )
        } else {
            figure + 
            axis_fnc(
                limits=limits,
                expand=expand,
                ...
            )
        }
    }
}

make_facet_formula <- function(
    figure,
    facet.group=NULL,
    facet.col=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.ncol=NULL,
    ...){
    if (!is.null(facet.col) & !is.null(facet.row)) {
        paste(
            paste(facet.row, collapse=' + '),
            paste(facet.col, collapse=' + '),
            sep=' ~ '
        ) %>%
        formula()
    } else if (!is.null(facet.row)) {
        formula(glue('{paste(facet.row, collapse=" + ")} ~ .'))
    } else if (!is.null(facet.col)) {
        formula(glue('~ {paste(facet.col, collapse=" + ")}'))
    } else if (!is.null(facet.group)) {
        formula(glue('~ {paste(facet.group, collapse=" + ")}'))
    } else {
        return(NULL)
        # stop(glue('Invalid faceting variable selection: {facet.group}. {facet.col}, {facet.row}'))
    }
}

add_faceting <- function(
    figure,
    scales='fixed',
    space='fixed',
    independent=FALSE,
    ragged='none',
    drop=TRUE,
    axes='margins',
    # margins='none',
    # solo_line=TRUE,
    # axes=FALSE,
    # trim_blank=TRUE,
    facet.group=NULL,
    facet.col=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.ncol=NULL,
    ...){
    faceting.formula <- 
        make_facet_formula(
            facet.group=facet.group,
            facet.col=facet.col,
            facet.row=facet.row,
            facet.nrow=facet.nrow,
            facet.ncol=facet.ncol
        )
    # Facet as specified
    if (!is.null(facet.group)) {
        figure <- 
            figure +
            # facet_wrap(
            facet_wrap2(
                vars(!!sym(facet.group)),
                nrow=facet.nrow,
                ncol=facet.ncol,
                scales=scales, 
                axes=independent,
                # drop=drop,
                # space=space,
                ...
            )
    } else {
        if (ragged == 'r') {
            figure <- 
                figure +
                facet_ragged_rows(
                    cols=vars(!!sym(facet.col)),
                    rows=vars(!!sym(facet.row)),
                    scales=scales,
                    axes=axes,
                    ...
                )
        } else if  (ragged == 'c') {
            figure <- 
                figure +
                facet_ragged_cols(
                    cols=vars(!!sym(facet.col)),
                    rows=vars(!!sym(facet.row)),
                    scales=scales,
                    axes=axes,
                    ...
                )
        } else {
            figure <- 
                figure +
                facet_nested(
                    faceting.formula,
                    scales=scales, 
                    space=space,
                    independent=independent,
                    drop=drop,
                    # solo_line=solo_line,
                    # axes=axes,
                    # trim_blank=trim_blank,
                    ...
                )
        }
    }
    figure
}

post_process_plot <- function(
    figure,
    plot.elements=NULL,
    theme.obj=NULL,
    scales='fixed',
    space='fixed',
    facet.ragged='none',
    independent=FALSE,
    drop=TRUE,
    # axes='margins',
    margins=FALSE,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.col=NULL,
    facet.ncol=NULL,
    facet.group=NULL,
    x.scale.mode='',
    x.log.base=10,
    x.axis.label.accuracy=0.1,
    x.n.breaks=NULL,
    x.limits=NULL,
    x.expand=c(0.00, 0.00, 0.00, 0.00),
    y.scale.mode='',
    y.log.base=10,
    y.axis.label.accuracy=0.01,
    y.n.breaks=NULL,
    y.limits=NULL,
    y.expand=c(0.00, 0.00, 0.00, 0.00),
    color.scale.mode='',
    color.log.base=10,
    color.axis.label.accuracy=0.1,
    color.n.breaks=NULL,
    color.limits=NULL,
    color.expand=c(0.00, 0.00, 0.00, 0.00),
    fill.scale.mode='',
    fill.log.base=10,
    fill.axis.label.accuracy=0.1,
    fill.n.breaks=NULL,
    fill.limits=NULL,
    fill.expand=c(0.00, 0.00, 0.00, 0.00),
    # linetype.scale.mode='',
    # linetype.log.base=10,
    # linetype.axis.label.accuracy=0.1,
    # linetype.n.breaks=NULL,
    # linetype.limits=NULL,
    # linetype.expand=c(0.00, 0.00, 0.00, 0.00),
    ...){
    figure %>% 
    add_faceting(
        scales=scales,
        space=space,
        independent=independent,
        ragged=facet.ragged,
        # axes=axes,
        drop=drop,
        # margins=margins,
        facet.row=facet.row,
        facet.col=facet.col,
        facet.group=facet.group,
        facet.nrow=facet.nrow,
        facet.ncol=facet.ncol
    ) %>% 
    # Set x-axis scaling (log, Mb, percent etc.)
    scale_axis(
        scale.axis='x',
        scale.mode=x.scale.mode,
        log.base=x.log.base,
        axis.label.accuracy=x.axis.label.accuracy,
        n.breaks=x.n.breaks,
        limits=x.limits,
        expand=x.expand
    ) %>%
    # Set y-axis scaling (log, Mb, percent etc.)
    scale_axis(
        scale.axis='y',
        scale.mode=y.scale.mode,
        log.base=y.log.base,
        axis.label.accuracy=y.axis.label.accuracy,
        n.breaks=y.n.breaks,
        limits=y.limits,
        expand=y.expand
    ) %>% 
    # Set color scaling (log, Mb, percent etc.)
    {
        if ('color' %in% names(figure@mapping)) {
            scale_axis(
                figure=.,
                scale.axis='color',
                scale.mode=color.scale.mode,
                log.base=color.log.base,
                axis.label.accuracy=color.axis.label.accuracy,
                n.breaks=color.n.breaks,
                limits=color.limits,
                expand=color.expand
            )
        } else {
            .
        }
    } %>% 
    # Set fill scaling (log, Mb, percent etc.)
    {
        if ('fill' %in% names(figure@mapping)) {
            scale_axis(
                figure=.,
                scale.axis='fill',
                scale.mode=fill.scale.mode,
                log.base=fill.log.base,
                axis.label.accuracy=fill.axis.label.accuracy,
                n.breaks=fill.n.breaks,
                limits=fill.limits,
                expand=fill.expand
            )
        } else {
            .
        }
    } %>% 
    # Add theme elements, as either an object or individual args
    {
        if (!is.null(theme.obj)) {
            . + theme.obj
        } else {
            . + make_ggtheme() + theme(...)
        } 
    } %>% 
    # Add extra elements
    {
        if (!is.null(plot.elements)) {
            purrr::reduce(
                .x=plot.elements,
                .f=function(x, y) { x + y },
                .init=.
            )
        } else {
            .
        }
    }
}

###################################################
# Modified quartabs render_tabset() function
###################################################
make_tabset_div <- function(
    level,
    pills=FALSE,
    tabset_width="default") {
    stopifnot(
    "`pills` must be a single logical value." =
        is.logical(pills) &&
            length(pills) == 1L &&
            !is.na(pills),
    "`tabset_width` must be one of 'default', 'fill', or 'justified'." =
        is.character(tabset_width) &&
            length(tabset_width) == 1L &&
            !is.na(tabset_width) &&
            tabset_width %in% c("default", "fill", "justified")
    )
    classes <- ".panel-tabset"
    if (isTRUE(pills)) {
    classes <- c(classes, ".panel-tabset-pills")
    }
    if (identical(tabset_width, "fill")) {
    classes <- c(classes, ".panel-fill")
    }
    if (identical(tabset_width, "justified")) {
    classes <- c(classes, ".panel-tabset-justified")
    }
    fence <- rep(':', level) %>% paste0(collapse='')
    paste0(fence, " {", paste(classes, collapse = " "), "}")
    # paste0("::: {", paste(classes, collapse = " "), "}")
}

make_fence <- function(level, bottom=FALSE){
    fence <- strrep(":", level)
    if (bottom) {
        cat(fence, "\n\n", sep="")
    } else {
        cat(fence, " {.panel-tabset}\n\n", sep="")
    }
}

make_heading <- function(label, level){
    stopifnot(
        "`level` must be an integer from 1 through 7." =
            length(level) == 1L &&
            !is.na(level) &&
            level == as.integer(level) &&
            level >= 1L &&
            level <= 7L
    )
    header <- 
        paste0(
            strrep("#", as.integer(level)),
            " ",
            as.character(label)
        )
    cat(header, '\n\n', sep="")
}

render_tabsets_recursively <- function(
    data,
    tabset_vars,
    output_var,
    hlevel,
    flevel,
    ...){
    current.variable <- tabset_vars[1]
    # recursive case, keep descending + printing titles until the bottom is reached
    if (length(tabset_vars) == 1) {
        make_fence(flevel)
        # pmap(
        #     .l=list(current.value=unique(data[[current.variable]])),
        data[[current.variable]] %>%
        { if (is.factor(.)) { levels(.) } else { unique(.) } } %>% 
        list() %>% 
        pmap(
            .l=.,
            .f=
                function(current.value) {
                    make_heading(current.value, hlevel)
                    data %>% 
                    filter(!!sym(current.variable) == current.value) %>% 
                    pull(output_var) %>%
                    {.[[1]]} %>% 
                    print()
                    cat('\n\n')
                }
        )
        make_fence(flevel, bottom=TRUE)
    } else {
        make_fence(flevel)
        # pmap(
        #     .l=list(current.value=unique(data[[current.variable]])),
        data[[current.variable]] %>%
        { if (is.factor(.)) { levels(.) } else { unique(.) } } %>% 
        list() %>% 
        pmap(
            .l=.,
            .f=
                function(current.value, data, current.variable, hlevel, ...){
                    make_heading(current.value, hlevel)
                    data %>% 
                    filter(!!sym(current.variable) == current.value) %>% 
                    render_tabsets_recursively(
                        hlevel=hlevel+1,
                        ...
                    )
                },
            current.variable=current.variable,
            data=data,
            tabset_vars=tabset_vars[2:length(tabset_vars)],
            output_var=output_var,
            flevel=flevel-1,
            hlevel=hlevel,
            ...
        )
        make_fence(flevel, bottom=TRUE)
    }
}

custom_render_tabset <- function(
    data,
    tabset_vars,
    output_var='figures',
    starting.hlevel=2,
    ...){
    # get nesting info
    current.variable <- tabset_vars[1]
    starting.fence.level <- length(tabset_vars) - 1 + 3
    # print quarto section headers 
    cat('\n\n')
    make_fence(starting.fence.level)
    # pmap(
    #     .l=.,
    data[[current.variable]] %>%
    { if (is.factor(.)) { levels(.) } else { unique(.) } } %>% 
    list() %>% 
    pmap(
        .l=.,
        .f=
            function(current.value, data, current.variable, hlevel, ...){
                make_heading(current.value, hlevel) 
                data %>% 
                filter(!!sym(current.variable) == current.value) %>% 
                render_tabsets_recursively(
                    hlevel=hlevel+1,
                    ...
                )
            },
        data=data,
        tabset_vars=tabset_vars[2:length(tabset_vars)],
        output_var=output_var,
        current.variable=current.variable,
        flevel=starting.fence.level-1,
        hlevel=starting.hlevel,
        ...
    )
    make_fence(starting.fence.level, bottom=TRUE)
    cat('\n\n')
}

###################################################
# Basic Plots
###################################################
plot_barplot <- function(
    plot.df,
    x.var=NULL,
    y.var=NULL,
    fill.var=NULL, 
    label.var=NULL,
    label.color='black',
    label.size=3,
    position='dodge',
    legend.cols=NA,
    ...){
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    # make it a boxplot 
    { 
        if (!is.na(legend.cols)) {
            . + 
            geom_col(position=position) +
            guides(fill=guide_legend(ncol=legend.cols))
        } else {
            . + 
            geom_col(position=position)
        }
    } %>% 
    # add text lables to bars (usually pcts)
    {
        if (!(is.null(label.var))) {
            . +
            geom_text(
                aes(label=.data[[label.var]]), 
                position=position_stack(vjust=0.5),
                # position=
                #     ifelse(
                #         position == 'dodge',
                #         position_dodge(),
                #         position_stack(vjust=0.5)
                #     ),
                color=label.color,
                size=label.size,
            )
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_boxplot <- function(
    plot.df,
    x.var=NULL,
    y.var=NULL,
    fill.var=NULL,
    outliers=FALSE,
    outlier.size=1,
    ...){
    # x.var='Comparison'; y.var='value'; fill.var='Edit.Numerator'; outliers=TRUE; outlier.size=1;
    {
        if (is.null(fill.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    fill=.data[[fill.var]]
                )
            )
        }
    } %>% 
    { . + geom_boxplot(outliers=TRUE, outlier.size=outlier.size) } %>% 
    # { . + geom_boxplot(outliers=FALSE) } %>% 
    # {
    #     if (outliers) {
    #         {.} + 
    #         geom_jitter(
    #             data=
    #                 .@data %>% 
    #                 filter(
    #                     !!sym(fill.var) < quantile(!!sym(fill.var), 0.25) - 1.5 * IQR(!!sym(fill.var)) | 
    #                     !!sym(fill.var) > quantile(!!sym(fill.var), 0.75) + 1.5 * IQR(!!sym(fill.var))
    #                 ) %>% 
    #                 ungroup(), 
    #             show.legend=FALSE,
    #             size=outlier.size
    #         )
    #     } else {
    #         .
    #     }
    # } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_density <- function(
    plot.df,
    x.var='',
    color.var=NULL, 
    alpha=0.2,
    binwidth=500,
    ...){
    # Set fill group if specified
    {
        # if (!is.null(fill.var)) {
        if (!is.null(color.var)) {
            ggplot(
                plot.df,
                aes(
                    fill=.data[[color.var]],
                    # color=.data[[color.var]],
                    x=.data[[x.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(x=.data[[x.var]])
            )
        }
    } %>% 
    # make it a boxplot 
    { . + geom_density(alpha=alpha) } %>% 
    # { . + geom_density_ridges(alpha=alpha) } %>% 
    # { . + geom_freqpoly(binwidth=binwidth, alpha=alpha) } %>% 
    # { . + geom_histogram(binwidth=binwidth, alpha=alpha) } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_ridges <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var=NULL,
    alpha=0.5,
    scales='fixed',
    ...){
    {
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
            )
        )
    } %>% 
    {
        if (!is.null(fill.var)) {
            . + aes(fill=.data[[fill.var]])
        } else {
            .
        }
    } %>% 
    {
        . + geom_density_ridges(alpha=alpha)
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_heatmap <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var='', 
    label.var=NULL,
    label.size=2,
    label.color='white',
    ...){
    {
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
                fill=.data[[fill.var]]
            )
        ) +
        geom_tile()
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        x.scale.mode='discrete',
        y.scale.mode='discrete',
        ...
    ) %>% 
    {
        if (!(is.null(label.var))) {
            . +
            geom_text(
                aes(label=.data[[label.var]]), 
                color=label.color,
                size=label.size,
            )
        } else {
            .
        }
    }
}

plot_jitter <- function(
    plot.df,
    x.var='',
    y.var='',
    color.var=NA, 
    shape.var=NA,
    regression_fnc=NULL,
    add_regression_SE=TRUE,
    alpha=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # Set fill group if specified
    # make it a scatter plot
    { 
        if (!is.na(color.var) & !is.na(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    color=.data[[color.var]],
                    shape=.data[[shape.var]]
                )
            )
        } else if (!is.na(color.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    color=.data[[color.var]]
                )
            )
        } else if (!is.na(shape.var)) {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]],
                    shape=.data[[shape.var]]
                )
            )
        } else {
            ggplot(
                plot.df,
                aes(
                    x=.data[[x.var]],
                    y=.data[[y.var]]
                )
            )
        }
    } %>% 
    {
        . +
        geom_jitter(
            alpha=alpha,
            size=size
        )
    } %>% 
    # add regression line(s) is specified
    {
        if (!is.null(regression_fnc)){
            . +
            geom_smooth(
                method=regression_fnc,
                se=add_regression_SE
            )
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_scatter_marginal <- function(
    plot.df,
    x.var='',
    y.var='',
    color.var=NA, 
    size.var=NA,
    regression_fnc=NULL,
    add_regression_SE=TRUE,
    alpha=0.5,
    size=0.5,
    margin.type='boxplot',
    show.outliers=FALSE,
    # which.margins='both',
    margin.size=0.3,
    scales='fixed',
    ...){
    is.color.categorical <- 
        case_when(
            is.na(color.var) ~ 'NA',
            is.factor(plot.df[[color.var]]) | is.character(plot.df[[color.var]]) ~ 'grouping',
            TRUE ~ 'continuous'
        )
    # specify ggplot varialbes based on whether color variable is null or categorical
    {
        if (is.color.categorical == 'grouping') {
            if (!is.na(size.var)) { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var, color=color.var, fill=color.var)) +
                geom_point(aes_string(size=size.var), alpha=alpha)
            } else { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var, color=color.var, fill=color.var)) +
                geom_point(size=size, alpha=alpha)
            }
        } else if (is.color.categorical == 'continuous') {
            if (!is.na(size.var)) { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var)) +
                geom_point(aes_string(color=color.var, size=size.var), alpha=alpha)
            } else { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var)) +
                geom_point(aes_string(color=color.var), size=size, alpha=alpha)
            }
        } else {
            if (!is.na(size.var)) { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var)) +
                geom_point(aes_string(size=size.var), alpha=alpha)
            } else { 
                ggplot(plot.df, aes_string(x=x.var, y=y.var)) +
                geom_point(size=size, alpha=alpha)
            }
        }
    } %>% 
    # add marginals to scatter plot
    {
        if (margin.type == 'density') {
            . +
            geom_xsidedensity(
                aes(y=after_stat(density)),
                alpha=alpha,
                position="stack"
            ) +
            geom_ysidedensity(
                aes(x=after_stat(density)),
                alpha=alpha,
                position="stack"
            ) +
            ggside(scales=scales)
        } else if (margin.type == 'boxplot') {
            . + 
            geom_xsideboxplot(
                orientation='x',
                scales=scales,
                outliers=show.outliers
            ) + 
            geom_ysideboxplot(
                orientation='y',
                scales=scales,
                outliers=show.outliers
            ) +
            ggside(scales=scales)
        } else if (margin.type == 'hist') {
            . + geom_ysidehistogram() + geom_xsidehistogram() + ggside(scales=scales)
        } else {
            .
        }
    } %>% 
    # add regression line(s) is specified
    {
        if (!is.null(regression_fnc)) {
            if (is.color.categorical == 'grouping') {
                . +
                geom_smooth(
                    aes_string(color=color.var),
                    method=regression_fnc,
                    se=add_regression_SE
                )
            } else { 
                . +
                geom_smooth(
                    method=regression_fnc,
                    se=add_regression_SE
                )
            }
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ggside.panel.scale=margin.size,
        ggside.panel.grid=element_blank(),
        ggside.panel.background=element_blank(),
        ...
    )
}

plot_lineplot <- function(
    plot.df,
    x.var='',
    y.var='',
    group.var='',
    color.var=NULL, 
    linetype.var=NULL,
    alpha=0.5,
    size=0.5,
    linewidth=1,
    scales='fixed',
    ...){
    # x.var='start'; x.scale.mode='mb'; y.var='nesting.lvl'; group.var='Sample.Group'; color.var='Genotype'; shape.var='Genotype'; facet.row='Edit';
    # make it a lineplot plot
    {
        ggplot(
            plot.df,
            aes(
                group=.data[[group.var]],
                x=.data[[x.var]],
                y=.data[[y.var]]
            )
        )
    } %>% 
    {
        if (!is.null(color.var)) {
            . + aes(color=.data[[color.var]])
            # . + aes(color=.data[[color.var]], fill=.data[[color.var]])
        } else {
            .
        }
    } %>% 
    {
        if (!is.null(linetype.var)) {
            . + 
            aes(
                linetype=.data[[linetype.var]],
                shape=.data[[linetype.var]]
            )
        } else {
            .
        }
    } %>% 
    { 
        . + 
        geom_point(size=size, alpha=alpha) + 
        geom_line(linewidth=linewidth) 
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_ribbonplot <- function(
    plot.df,
    x.var='',
    y.var='',
    ymin.var='',
    ymax.var='',
    color.var=NULL, 
    alpha=0.5,
    size=0.5,
    linewidth=1,
    # scales='fixed',
    ...){
    # make it a lineplot plot
    {
        ggplot(
            plot.df,
            aes(
                # group=.data[[group.var]],
                x=.data[[x.var]],
                y=.data[[y.var]],
                ymin=.data[[ymin.var]],
                ymax=.data[[ymax.var]]
            )
        )
    } %>% 
    {
        if (!is.null(color.var)) {
            . + 
            geom_ribbon(aes(fill=.data[[color.var]]), alpha=alpha) +
            geom_point(aes(color=.data[[color.var]]), size=size) +
            geom_line(aes(color=.data[[color.var]]), linewidth=linewidth) 
        } else {
            . + 
            geom_ribbon(alpha=alpha) +
            geom_point(size=size) +
            geom_line(linewidth=linewidth) 
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        # scales=scales,
        ...
    )
}

