## @knitr report_initialization

library(showtext)

# Check the current search path for fonts
font_paths()
#> [1] "C:\\Windows\\Fonts"

# List available font files in the search path
f <- font_files()

font_add("CMU Classical Serif", "cmunci.ttf")
font_add("CMU Serif Upright Italic", "cmunui.ttf")
font_add("CMU Serif", "cmunrm.ttf")

font_add("Merriweather Regular", "Merriweather-Regular.ttf")
font_add("Merriweather Light", "Merriweather Light.ttf")

font_families()

## automatically use showtext for new devices
showtext_auto()


p.colors <-
  brewer.pal(name = "YlOrBr", n = 9)[c(5, 7, 9)]

palette_nar <- brewer.pal(name = "YlOrBr", n = 9)[c(5, 7, 9)]
palette_chiar <- brewer.pal(name = "GnBu", n = 9)[c(5, 7, 9)]
palette_binar <- brewer.pal(name = "YlGn", n = 9)[c(5, 7, 9)]
palette_podar <- brewer.pal(name = "BuPu", n = 9)[c(5, 7, 9)]

# Pair plots --------------------------------------------------------------

## @knitr harvest_cleanup

harvest_cleanup <- function(harv) {

  d <- harv %>%
    filter(standardization == "unstd",
           param.name == "X.WITH.PHI") %>%
    select(type,
           N,
           T,
           l2.dist,
           Model,
           est,
           sig,
           lower_2.5ci,
           upper_2.5ci) %>%
    group_by(N, T, type, Model, l2.dist) %>%
    mutate(sign.X.sig = as.factor(sign(est) * sig)) %>%
    mutate(
      mean.est = mean(est),
      n.datasets = n(),
      nonconverged.percent = round(100 * (1000 - n()) / 1000, 2)
    ) %>%
    group_by(sign.X.sig,
             .add = TRUE) %>%
    mutate(error.percents = round(100 * n() / n.datasets, 2)) %>%
    arrange(N, T,
            .by_group = TRUE) %>%
    mutate(ordering = est) %>%
    arrange(ordering) %>%
    ungroup() %>%
    group_by(type,
             N,
             T,
             l2.dist,
             Model) %>%
    mutate(ord = order(ordering) - n() / 2) %>%
    mutate(NN = as.factor(paste("N =", N)),
           TT = as.factor(paste("T =", T)))


  levels(d$sign.X.sig) <- c("Significant negative estimates",
                            "Non-significant estimates",
                            "Significant positive estimates")


  d.abridged <- harv %>%
    na.omit() %>%
    mutate(std_parname = paste(standardization, param.name)) %>%
    filter(std_parname %in% c("stdyx X.WITH.PHI",
                              "unstd X.WITH.PHI",
                              "unstd Means.PHI")) %>%
    select(-standardization) %>%
    mutate(
      parameter = case_when(
        std_parname  == "stdyx X.WITH.PHI" ~ "Correlation",
        std_parname  == "unstd X.WITH.PHI" ~ "Covariance",
        std_parname  == "unstd Means.PHI" ~ "Fixed.Phi",
      ),
      .after = Rep
    ) %>%
    select(uSeed:sig) %>%
    mutate(
      true.value = case_when(
        parameter == "Correlation" ~ 0,
        parameter == "Covariance" ~ 0,
        parameter == "Fixed.Phi" ~ 0.4
      ),
      .after = sig
    )


  d_summary <- d.abridged %>%
    select(type:nonconverged.percent) %>%
    select(-Rep, -est) %>%
    distinct() %>%
    mutate(`Type-1 error` = percent,
           .after = Variance)


  dd <- d_summary[with(d_summary,
                       order(type,
                             l2.X.Model,
                             parameter,
                             N,
                             T)), ] %>%
    ungroup() %>%
    complete(Resid, `Model name`, N, T,
             fill = list(value = 0))


  d_important <- dd %>%
    mutate(
      `Model name` = factor(
        paste0(Model, "(1)"),
        levels = c("NAR(1)",
                   "Chi2AR(1)",
                   "BinAR(1)",
                   "PoDAR(1)")
      ),
      `Means distribution` = paste0(l2.dist, "-distributed means"),
      Resid = case_when(
        type == "resid.fixed" ~ "`Fixed Residual Variance`",
        type == "resid.random" ~ "`Random Residual Variance`"
      ),
      .after = l2.X.Model
    )

  levels(d_important$`Model name`) <- c(
    `NAR(1)` = "AR(1)",
    `Chi2AR(1)` = TeX("$\\chi^2$AR(1)"),
    `BinAR(1)` = "BinAR(1)",
    `PoDAR(1)` = "PoDAR(1)"
  )




}


## @knitr plot_histograms

plot_histograms <- function(d,
                            binwidth = 0.5,
                            nrow.facet = 20,
                            p.colors = c("#FE9929",
                                         "#CC4C02",
                                         "#662506")) {
  fill_col <- p.colors[2]
  x.range <- d$x %>% range()

  d %>%
    group_by(subject) %>%
    mutate(mm = mean(x)) %>%
    arrange(desc(mm)) %>%
    ggplot(aes(x = x,
               y = ..ndensity..)) +
    geom_histogram(binwidth = binwidth,
                   fill = fill_col,
                   center = 0) +
    facet_wrap(~ reorder(subject, mm),
               nrow = nrow.facet) +
    geom_hline(yintercept = 0,
               size = 0.1) +
    ggtitle(NULL) +
    theme_tufte() +
    theme(
      strip.background = element_blank(),
      aspect.ratio = 1,
      strip.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

}

## @knitr plot_pairplots

plot_pairplots <- function(d,
                           p.colors = c("#FE9929",
                                        "#CC4C02",
                                        "#662506"),
                           rel_title = 2,
                           rel_lines = 0.5,
                           rel_dots = 3,
                           bins = 30) {
  color_hist <- p.colors[2]
  color_scatter <- p.colors[1]

  d <- d %>%
    group_by(subject) %>%
    summarise(
      Mean = mean(x,
                  na.rm = TRUE),
      Variance = var(x,
                     na.rm = TRUE),
      Skewness = moments::skewness(x,
                                   na.rm = TRUE)
    ) %>%
    na.omit()

  cors <- d %>%
    select(Mean, Variance, Skewness) %>%
    correlation::correlation() %>%
    as.data.frame() %>%
    mutate(
      value = round(r, 3),
      significance = cut(
        p,
        breaks = c(0, 0.001, 0.01, 0.05, 1),
        include.lowest = T,
        labels = c('***', '**', '*', '')
      ),
      string = paste0(value, significance)
    ) %>%
    select(-r:-n_Obs)

  margin.thing <- 5

  ##  Correlations

  p_cor_Mean.Variance <-
    ggplot() +
    geom_text(aes(
      x = mean(range(d$Variance)),
      y = mean(range(d$Mean)),
      label = paste(
        "Corr:",
        cors %>%
          filter(Parameter1 == "Mean",
                 Parameter2 == "Variance") %>%
          pull(string)
      )
    )) +
    theme_tufte() +
    theme(
      text = element_text(size = rel(rel_title)),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  p_cor_Mean.Skewness <-
    ggplot() +
    geom_text(aes(
      x = mean(range(d$Skewness)),
      y = mean(range(d$Mean)),
      label = paste(
        "Corr:",
        cors %>%
          filter(Parameter1 == "Mean",
                 Parameter2 == "Skewness") %>%
          pull(string)
      )
    )) +
    theme_tufte() +
    theme(
      text = element_text(size = rel(rel_title)),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  p_cor_Variance.Skewness <-
    ggplot() +
    geom_text(aes(
      x = mean(range(d$Skewness)),
      y = mean(range(d$Variance)),
      label = paste(
        "Corr:",
        cors %>%
          filter(Parameter1 == "Variance",
                 Parameter2 == "Skewness") %>%
          pull(string)
      )
    )) +
    theme_tufte() +
    theme(
      text = element_text(size = rel(rel_title)),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  ## Histograms

  p_mean <-
    d %>%
    ggplot(aes(x = Mean)) +
    geom_histogram(fill = color_hist,
                   bins = bins) +
    ggtitle("Mean") +
    theme_tufte() +
    theme(
      plot.title = element_text(
        size = rel(rel_title),
        hjust = 0.5,
        margin = margin(t = margin.thing,
                        b = -margin.thing)
      ),
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )

  p_variance <-
    d %>%
    ggplot(aes(x = Variance)) +
    geom_histogram(fill = color_hist,
                   bins = bins) +
    ggtitle("Variance") +
    theme_tufte() +
    theme(
      plot.title = element_text(
        size = rel(rel_title),
        hjust = 0.5,
        margin = margin(t = margin.thing,
                        b = -margin.thing)
      ),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    )

  p_skewness <-
    d %>%
    ggplot(aes(x = Skewness)) +
    geom_histogram(fill = color_hist,
                   bins = bins) +
    ggtitle("Skewness") +
    geom_vline(
      xintercept = 0,
      alpha = 0.5,
      color = "gray45",
      linetype = "solid",
      size = rel(rel_lines)
    ) +
    geom_vline(
      xintercept = c(-0.5, 0.5),
      color = "gray45",
      linetype = "dotted",
      size = rel(rel_lines)
    ) +
    geom_vline(
      xintercept = c(-1, 1),
      color = "gray45",
      linetype = "dashed",
      size = rel(rel_lines)
    ) +
    theme_tufte() +
    theme(
      plot.title = element_text(
        size = rel(rel_title),
        hjust = 0.5,
        margin = margin(t = margin.thing,
                        b = -margin.thing)
      ),
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )

  ## scatter plots

  p_scatter_Mean.Variance <- d %>%
    ggplot(aes(x = Mean,
               y = Variance)) +
    geom_point(
      color = color_scatter,
      shape = "bullet",
      alpha = 0.85,
      size = rel(rel_dots)
    ) +
    theme_tufte() +
    theme(
      axis.title = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )



  p_scatter_Mean.Skewness <- d %>%
    ggplot(aes(x = Mean,
               y = Skewness)) +
    geom_hline(
      yintercept = 0,
      alpha = 0.5,
      color = "gray45",
      linetype = "solid",
      size = rel(rel_lines)
    ) +
    geom_hline(
      yintercept = c(-0.5, 0.5),
      color = "gray45",
      linetype = "dotted",
      size = rel(rel_lines)
    ) +
    geom_hline(
      yintercept = c(-1, 1),
      color = "gray45",
      linetype = "dashed",
      size = rel(rel_lines)
    ) +
    geom_point(
      color = color_scatter,
      shape = "bullet",
      alpha = 0.7,
      size = rel(rel_dots)
    ) +
    theme_tufte() +
    theme(axis.title = element_blank())



  p_scatter_Variance.Skewness <- d %>%
    ggplot(aes(x = Variance,
               y = Skewness)) +
    geom_hline(
      yintercept = 0,
      alpha = 0.5,
      color = "gray45",
      linetype = "solid",
      size = rel(rel_lines)
    ) +
    geom_hline(
      yintercept = c(-0.5, 0.5),
      color = "gray45",
      linetype = "dotted",
      size = rel(rel_lines)
    ) +
    geom_hline(
      yintercept = c(-1, 1),
      color = "gray45",
      linetype = "dashed",
      size = rel(rel_lines)
    ) +
    geom_point(
      color = color_scatter,
      shape = "bullet",
      alpha = 0.7,
      size = rel(rel_dots)
    ) +
    theme_tufte() +
    theme(
      axis.title = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank()
    )

  ## Putting plots together

  p_pairplots <-
    (p_mean + p_cor_Mean.Variance + p_cor_Mean.Skewness) /
    (p_scatter_Mean.Variance + p_variance + p_cor_Variance.Skewness) /
    (p_scatter_Mean.Skewness + p_scatter_Variance.Skewness + p_skewness) &
    # plot_layout(guides = "collect") &
    plot_annotation(theme = theme(
      plot.title = element_text(family = "Merriweather Regular"),
      axis.text = element_text(size = rel(rel_title) / 2)
    ))

  ## Returning the 3x3 pairplots
  p_pairplots

}


## @knitr plot_dataset.profile

plot_dataset.profile <- function(sim.object,
                                 p.colors = NULL,
                                 title_l2.dist = FALSE) {
  d_sim <- sim.object$sim.Dataset

  l2.distr <- ifelse(sim.object$l2.dist == "Chi2",
                     "$\\chi^2$",
                     "Gaussian")
  model <- sim.object$Model
  title <- model

  binwidth <- 0.5
  if (model == "NAR") {
    title <- "AR"
    palette_color <- palette_nar
  }
  if (model == "Chi2AR") {
    title <- "$\\chi^2$AR"
    palette_color <- palette_chiar
  }
  if (model == "BinAR")
    palette_color <- palette_binar
  if (model == "PoDAR")
    palette_color <- palette_podar

  if (is.null(p.colors))
    p.colors <- palette_color

  title <- paste0(title,
                  "(1)",
                  " $\\phantom{\\chi^2}$")

  if (title_l2.dist)
    title <- paste0(l2.distr,
                    "-distributed means",
                    " $\\phantom{\\chi^2}$")


  p_out <-  plot_histograms(
    d = d_sim,
    binwidth = binwidth,
    nrow.facet = 10,
    p.colors = p.colors
  ) +
    plot_spacer() +
    plot_pairplots(d = d_sim,
                   p.colors = p.colors) +
    plot_layout(widths = c(1, 0.05, 1)) +
    plot_annotation(
      title = TeX(title),
      theme = theme(
        plot.title = element_text(size = rel(2.5),
                                  family = "CMU Serif"),
        plot.subtitle = element_blank()
      )
    )

  p_out %>% wrap_elements()

}


## @knitr save_dataset_profile

save_dataset_profile <-
  function(upper,
           lower,
           file.name = "meh",
           title = NULL) {
    p_profiles <-
      (
        plot_dataset.profile(upper, title_l2.dist = TRUE) /
          plot_spacer() /
          plot_dataset.profile(lower, title_l2.dist = TRUE)
      ) +
      plot_layout(heights = c(1, 0.02, 1)) +
      plot_annotation(# title = ifelse(is.null(title),
        #                NULL,
        #                TeX(paste(title, "model $\\phantom{\\chi^2}$"))
        #                ),
        theme = theme(
          plot.title = element_text(
            size = rel(4),
            hjust = 0.5,
            family = "CMU Serif"
          ),
          plot.subtitle = element_blank()
        ))

    ggsave(
      paste0(file.name,
             ".pdf"),
      p_profiles,
      width = 2 * 15,
      height = 2.02 * 1.1 * 15,
      units = "cm"
    )

  }

levels(d$sign.X.sig) <- c(
  "Significant negative estimates",
  "Non-significant estimates",
  "Significant positive estimates"
)


## @knitr plot_Model.x.Resid

plot_Model.x.Resid <- function(d_important,
                               which.parameter,
                               which.measure = "Positive error",
                               Means.dist = c("Gaussian", "Chi2"),
                               line.width = 1,
                               font.scale = 5) {
  ddd <- d_important %>%
    filter(parameter == which.parameter)

  if (!grepl("error",
             tolower(which.measure),
             fixed = TRUE)) {
    ddd$value <- ddd %>% pull(which.measure)
    ddd <- ddd %>% filter(sign.X.sig == "Zero")
    y.axis <- which.measure
  }

  if (grepl("positive",
            tolower(which.measure),
            fixed = TRUE)) {
    which.measure <- "Type-1 error"
    ddd$value <- ddd %>% pull(which.measure)
    ddd <- ddd %>% filter(sign.X.sig == "Positive")
    y.axis <- "Positive error"
  }

  if (grepl("negative",
            tolower(which.measure),
            fixed = TRUE)) {
    which.measure <- "Type-1 error"
    ddd$value <- ddd %>% pull(which.measure)
    ddd <- ddd %>% filter(sign.X.sig == "Negative")
    y.axis <- "Negative error"
  }

  if (grepl("total",
            tolower(which.measure),
            fixed = TRUE)) {
    which.measure <- "Type-1 error"
    ddd$value <- 100 - (ddd %>% pull(which.measure))
    ddd <- ddd %>% filter(sign.X.sig == "Zero")
    y.axis <- "Total error"
  }


  y.range <- ddd$value %>% range()
  y.range[1] <- min(-0.001, y.range[1])
  y.range[2] <- max(0.001, y.range[2])

  ddd <- ddd  %>%
    filter(l2.dist == Means.dist)

  ddd$N <- as.factor(ddd$N)
  ddd$T <- as.factor(ddd$T)


  title <- ifelse(
    Means.dist == "Gaussian",
    TeX("Gaussian-distributed means $\\phantom{\\chi^2}$"),
    TeX("$\\chi^2$-distributed means")
  )


  output.plot <- ddd %>%
    ggplot() +
    aes(
      x = T,
      y = value,
      group = N,
      color = N
    ) +
    ## make solid y=0 axis line
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      color = "black",
      alpha = 1
    ) +
    geom_line(size = rel(1),
              alpha = 0.8,
              lineend = "round") +
    geom_point(size = rel(1.5), alpha = 1) +
    scale_color_manual(values = brewer.pal(name = "PuRd", n = 9)[c(4, 6, 9)]) +
    facet_grid(
      rows = vars(`Model name`),
      cols = vars(Resid),
      labeller = label_parsed
    ) +
    theme_light() +
    theme_pubclean() +
    scale_y_continuous(limits = y.range) +
    ylab(y.axis) +
    ggtitle(title) +
    theme(
      plot.title = element_text(
        family = "CMU Serif"),
      panel.spacing = unit(0.7, "lines"),
      legend.position = "bottom",
      legend.key = element_rect(colour = NA, fill = NA),
      legend.key.width = unit(10, "mm"),
      text = element_text(size = 10,
                          family = "CMU Serif")
    )

  if (which.measure == "Type-1 error" |
      which.measure == "Relative efficiency") {
    output.plot <- ddd %>%
      ggplot() +
      aes(
        x = T,
        y = value,
        group = N,
        color = N
      ) +
      ## make solid y=0 axis line
      geom_hline(
        yintercept = 0,
        linetype = "solid",
        color = "black",
        alpha = 1
      ) +
      ## manually add needed grids
      geom_hline(
        yintercept = c(10, 20, 40, 60),
        linetype = "dotted",
        color = "#BEBEBE",
        alpha = 1
      ) +
      geom_line(size = rel(1),
                alpha = 0.8,
                lineend = "round") +
      geom_point(size = rel(1.5), alpha = 1) +
      scale_color_manual(values = brewer.pal(name = "PuRd", n = 9)[c(4, 6, 9)]) +
      facet_grid(
        rows = vars(`Model name`),
        cols = vars(Resid),
        labeller = label_parsed
      ) +
      theme_light() +
      theme_pubclean() +
      scale_y_continuous(breaks = c(2.5, 10, 20, 40, 60, 80, 90, 100),
                         # looked up myself to assure  + and - ranges are equal
                         limits = c(0, 61)) +
      ylab(y.axis) +
      ggtitle(title) +
      theme(
        plot.title = element_text(#size = rel(1.5),# 4*font.scale,
          family = "CMU Serif"),
        panel.spacing = unit(0.7, "lines"),
        legend.position = "bottom",
        legend.key = element_rect(colour = NA, fill = NA),
        legend.key.width = unit(10, "mm"),
        # legend.key.height = unit(10, "mm"),
        text = element_text(size = 10,
                            family = "CMU Serif"),
        ## remove grid lines
        panel.grid.major.y = element_blank()
      )
    # add thresholds
    # geom_hline(yintercept = 2.5,
    #            linetype = "dotted",
    #            alpha = 0.7) +
    if (y.axis == "Total error")
      output.plot <- output.plot +
        geom_hline(yintercept = 5,
                   linetype = "dashed",
                   alpha = 0.7)

    if (y.axis != "Total error")
      output.plot <- output.plot +
        geom_hline(yintercept = 2.5,
                   linetype = "dashed",
                   alpha = 0.7)

  }

  output.plot <- output.plot +
    geom_hline(yintercept = 0,
               alpha = 0.6)

  return(output.plot)

}



## @knitr plot_quadrants

plot_quadrants <- function(d_important,
                           which.parameter = "Correlation",
                           upper.measure = "Bias",
                           lower.measure = "RMSE",
                           line.width = 2,
                           font.scale = 8) {
  p.upper.left <- plot_Model.x.Resid(d_important,
                                     which.parameter,
                                     upper.measure,
                                     "Gaussian",
                                     line.width,
                                     font.scale) + xlab(NULL)

  p.upper.right <- plot_Model.x.Resid(d_important,
                                      which.parameter,
                                      upper.measure,
                                      "Chi2",
                                      line.width,
                                      font.scale) + xlab(NULL) + ylab(NULL)


  p.lower.left <- plot_Model.x.Resid(d_important,
                                     which.parameter,
                                     lower.measure,
                                     "Gaussian",
                                     line.width,
                                     font.scale)  + ggtitle(NULL)

  p.lower.right <- plot_Model.x.Resid(d_important,
                                      which.parameter,
                                      lower.measure,
                                      "Chi2",
                                      line.width,
                                      font.scale)  + ggtitle(NULL) + ylab(NULL)


  p.patchwork <-
    (p.upper.left | p.upper.right) / (p.lower.left | p.lower.right)

  title <- paste(upper.measure,
                 "and",
                 lower.measure,
                 "in the estimated",
                 tolower(which.parameter))

  if (which.parameter == "Fixed.Phi")
    title <- "Level-2 $\\phi$"

  if (grepl("error", tolower(upper.measure), fixed = TRUE)) {
    title <- paste("Type-I error rates in the estimated",
                   tolower(which.parameter))
  }


  p.final <- p.patchwork +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(plot.title =
                      element_text(
                        size = 20,
                        family = "CMU Serif",
                        hjust = 0.5
                      ))) &
    theme(legend.position = "bottom")

  return(p.final)

}


## @knitr plot_caterpillar

plot_caterpillar <- function(d,
                             l2.dist_ = "Chi2",
                             Model_ = "PoDAR",
                             analysis.type = "resid.fixed",
                             parameter = "covariance",
                             legend.key.width = 5,
                             legend.line.width = 10) {
  title <- ifelse(
    analysis.type == "resid.fixed",
    "fixed residual variance",
    "random residual variance"
  ) %>%
    paste("Modeled with", .)

  dd <- d %>%
    filter(l2.dist == l2.dist_,
           Model == Model_)
  y.range <- c(min(dd$lower_2.5ci), max(dd$upper_2.5ci))

  dd %>%
    filter(type == analysis.type) %>%
    # sample_n(100) %>%
    ggplot() +
    aes(x = ord,
        color = sign.X.sig) +
    geom_segment(
      aes(
        x = ord,
        xend = ord,
        y = upper_2.5ci,
        yend = lower_2.5ci
      ),
      alpha = 0.8 + 0.2,
      size = rel(0.06)
    ) +
    geom_segment(
      aes(
        x = ord,
        xend = ord,
        y = est + min(0.01, 0.05 * abs(upper_2.5ci - lower_2.5ci)),
        yend = est - min(0.01, 0.05 * abs(upper_2.5ci - lower_2.5ci))
      ),
      alpha = 0.8 + 0.2,
      color = "white",
      size = rel(0.06)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "solid",
      size = rel(0.15),
      color = "black",
      alpha = 0.8
    ) +
    theme_pubclean() +
    guides(colour = guide_legend(override.aes = list(
      linewidth = rel(1.5),
      alpha = 1
    ))) +
    scale_color_manual(values =
                         c("#F67E4BFF",
                           "#98CAE1FF",
                           "#364B9AFF")) +
    scale_y_continuous(limits = y.range) +
    facet_grid(rows = vars(NN),
               # rows = vars(rev(TT)),
               cols = vars(TT)) +
    ggtitle(title) +
    theme(
      legend.background = element_rect(colour = NA, fill = NA),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      # "Score",
      axis.text.x = element_blank(),
      legend.key.width = unit(legend.key.width * 2, "mm"),
      axis.text.y = element_text(size = 10),
      text = element_text(family = "CMU Serif",
                          size = 12)
    ) +
    xlab(NULL) +
    ylab(paste("Estimated", parameter))

}
