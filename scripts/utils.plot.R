###################################################
# Dependencies
###################################################
suppressPackageStartupMessages({
    library(ggplot2)
    library(ggpubr)
    library(ggh4x)
    # library(ggridges)
    library(GGally)
    library(scales)
    # library(ggpointdensity)
    # library(viridis)
    library(ComplexUpset)
    library(quarto)
    library(quartabs)
})

###################################################
# Transform data for plotting
###################################################
calc_pct <- function(
    count.df,
    cols_exclude=c(),
    col_pct=NULL){
    # count.df=n.TADCompare.df; cols_exclude=c('chr', 'DifferenceType', 'Enriched.Condition'); col_pct=NULL
    # calculate relative frequency from count data
    count.df %>% 
    {
        if ('n' %in% colnames(.)){
            group_by(., across(-c(cols_exclude, 'n'))) %>% 
            summarize(n=sum(n))
        } else {
            group_by(., across(-c(cols_exclude))) %>% 
            count()
        }
    } %>% 
    {
        if (!is.null(col_pct)) {
            ungroup(.) %>%
            group_by(across(-c(col_pct, 'n')))
        } else {
            .
        }
    } %>% 
    mutate(total=sum(n)) %>% 
    ungroup() %>% 
    mutate(pct=n / total) %>% 
    mutate(n.label=glue('n = {n}')) %>% 
    mutate(pct.label=glue('{round(100 * n / total, digits=1)}%')) %>% 
    mutate(n.and.pct.label=glue('{pct.label}\n({n.label})'))
}

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

add_faceting <- function(
    figure,
    scales='fixed',
    space='fixed',
    independent=FALSE,
    # solo_line=TRUE,
    # axes=FALSE,
    # trim_blank=TRUE,
    facet.group=NULL,
    facet.col=NULL,
    facet.row=NULL,
    facet.nrow=NULL,
    facet.ncol=NULL,
    ...){
    # Facet as specified
    if (!is.null(facet.col) & !is.null(facet.row)) {
        figure <- 
            figure +
            facet_nested(
                paste(
                    paste(facet.row, collapse=' + '),
                    paste(facet.col, collapse=' + '),
                    sep=' ~ '
                ) %>%
                formula(),
                scales=scales, 
                space=space,
                independent=independent,
                # solo_line=solo_line,
                # axes=axes,
                # trim_blank=trim_blank,
                ...
            )
    } else if (!is.null(facet.row)) {
        figure <- 
            figure +
            facet_nested(
                formula(glue('{paste(facet.row, collapse=" + ")} ~ .')),
                scales=scales, 
                space=space,
                independent=independent,
                # solo_line=solo_line,
                # axes=axes,
                # trim_blank=trim_blank,
                ...
            )
    } else if (!is.null(facet.col)) {
        figure <- 
            figure +
            facet_nested(
                formula(glue('~ {paste(facet.col, collapse=" + ")}')),
                scales=scales, 
                space=space,
                independent=independent,
                # solo_line=solo_line,
                # trim_blank=trim_blank,
                # axes=axes,
                ...
            )
    } else if (!is.null(facet.group)) {
        figure <- 
            figure +
            facet_wrap2(
                vars(!!sym(facet.group)),
                nrow=facet.nrow,
                ncol=facet.ncol
                # ...
            )
    }
    figure
}

post_process_plot <- function(
    figure,
    plot.elements=NULL,
    theme.obj=NULL,
    scales='fixed',
    space='fixed',
    independent=FALSE,
    drop=TRUE,
    axes='margins',
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
    y.axis.label.accuracy=0.1,
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
        axes=axes,
        drop=drop,
        margins=margins,
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
    # Set linetype scaling (log, Mb, percent etc.)
    # {
    #     if ('linetype' %in% names(figure@mapping)) {
    #         scale_axis(
    #             figure=.,
    #             scale.axis='linetype',
    #             scale.mode=linetype.scale.mode,
    #             log.base=linetype.log.base,
    #             axis.label.accuracy=linetype.axis.label.accuracy,
    #             n.breaks=linetype.n.breaks,
    #             limits=linetype.limits,
    #             expand=linetype.expand
    #         )
    #     } else {
    #         .
    #     }
    # } %>% 
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

plot_histogram <- function(
    plot.df,
    x.var='',
    fill.var=NULL, 
    alpha=0.7,
    binwidth=0.05,
    position='stack',
    ...){
    { 
        ggplot(
            plot.df,
            aes(x=.data[[x.var]])
        ) 
    } %>%
    {
        if (!is.null(fill.var)) {
            . + aes(fill=.data[[fill.var]])
            # . + aes(color=.data[[fill.var]])
        } else {
            .
        }
    } %>% 
    {
        . +
        # geom_freqpoly(
        geom_histogram(
            position=position,
            alpha=alpha,
            binwidth=binwidth
        ) 
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
}

plot_violin <- function(
    plot.df,
    x.var='',
    y.var='',
    fill.var=NULL, 
    plot_pts=FALSE,
    jitter.size=0.1,
    quantile.color=NULL,
    draw.quantiles=0L,
    quantile.linewidth=NULL,
    position='dodge',
    adjust=0.5,
    ...){
    # make it a violin plot
    { 
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]],
                fill=.data[[fill.var]],
            )
        )
    } %>%
    {
        .+
        geom_violin(
            quantile.color=quantile.color,
            quantile.linetype=draw.quantiles,
            quantile.linewidth=quantile.linewidth,
            position=position,
            adjust=adjust
        )
    } %>% 
    # plot individual points if specified
    {
        if (plot_pts){
            . +
            geom_jitter(size=jitter.size)
            # geom_jitter(aes(size=jitter.size))
        } else {
            .
        }
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(...)
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
    { . + geom_line(linewidth=linewidth) } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        ...
    )
}

plot_pointdensity <- function(
    plot.df,
    x.var='',
    y.var='',
    adjust=0.5,
    size=0.5,
    scales='fixed',
    ...){
    # make it a scatter plot
    { 
        ggplot(
            plot.df,
            aes(
                x=.data[[x.var]],
                y=.data[[y.var]]
            )
        ) +
        geom_pointdensity(
            adjust=adjust,
            size=size
        ) +
        scale_color_viridis()
    } %>% 
    # Handle faceting + scaling + theme options
    post_process_plot(
        scales=scales,
        color.scale.mode='',
        fill.scale.mode='',
        ...
    )
}

plot_ggpairs <- function(
    freq.df,
    group.colname='Sample.Group',
    value.colname='value',
    color.colname=NULL,
    pt.alpha=0.6,
    text.size=8,
    linewidth=0.15,
    ...){
    plot.df <- 
        freq.df %>% 
        pivot_wider(
            names_from=!!sym(group.colname),
            names_prefix='group_',
            values_from=!!sym(value.colname)
        ) %>%  
        unnest(starts_with('group_'))
    cols.to.keep <- 
        plot.df %>% 
        summarize(across(everything(), ~ sum(!is.na(.x)) > 0)) %>%
        pivot_longer(everything(), names_to='colname', values_to='is.not.all.NAs') %>% 
        mutate(is.group.col=str_detect(colname, '^group_')) %>%
        filter((is.group.col & is.not.all.NAs) | !is.group.col) %>% 
        pull(colname) %>% 
        c(., color.colname) %>%
        unique()
    plot.df <- 
        plot.df %>% 
        select(cols.to.keep)
    cols <- 
        grep('^group_', colnames(plot.df), value=TRUE)
    cols.names <- 
        str_remove(cols, '^group_')
    mapping <- 
        if (!is.null(color.colname)) {
            aes(color=.data[[color.colname]], alpha=pt.alpha)
        } else {
            aes(alpha=pt.alpha)
        }
    # print(cols.to.keep)
    # return(plot.df)
    ggpairs(
        plot.df,
        mapping=mapping,
        # columns=which(colnames(plot.df) %in% cols),
        columns=cols,
        columnLabels=cols.names
    ) +
    geom_abline(intercept=0, slope= 1, linewidth=linewidth, color='black', linetype='dashed') +
    geom_abline(intercept=0, slope=-1, linewidth=linewidth, color='black', linetype='dashed') +
    geom_vline(xintercept=0,           linewidth=linewidth, color='black', linetype='solid') +
    geom_hline(yintercept=0,           linewidth=linewidth, color='black', linetype='solid') +
    theme(
        axis.text.x=element_text(size=text.size),
        axis.text.y=element_text(size=text.size),
        strip.text=element_text(size=text.size+1)
    ) + 
    make_ggtheme()
}

