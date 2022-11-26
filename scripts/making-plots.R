## @knitr results_read_harvest

d_abridged <- list.files(here::here(dir_harvest),
                          pattern = "harvest-abridged_") %>%
  sort(decreasing = TRUE) %>%
  paste(here::here(dir_harvest),
        .,
        sep = "/") %>%
  read_rds()

d_important <- list.files(here::here(dir_harvest),
                          pattern = "harvest-important_") %>%
  sort(decreasing = TRUE) %>%
  paste(here::here(dir_harvest),
        .,
        sep = "/") %>%
  read_rds()


## @knitr results_plots_measures_single

for (parameter in c("Correlation", "Covariance")) {
  for (outcome in c("Bias",
                    "Variance",
                    "RMSE",
                    "MAE",
                    "Positive error",
                    "Negative error",
                    "Total error")) {


    p.error.left <- plot_Model.x.Resid(d_important,
                                       parameter,
                                       outcome,
                                       "Gaussian",
                                       line.width,
                                       font.scale)


    p.error.right <- plot_Model.x.Resid(d_important,
                                        parameter,
                                        outcome,
                                        "Chi2",
                                        line.width,
                                        font.scale) + ylab(NULL)


    p.patchwork <- (p.error.left | p.error.right)


    title <- paste(outcome,
                   "in the estimated",
                   tolower(parameter)
    )

    p.final <- p.patchwork +
      plot_layout(guides = "collect") +
      plot_annotation(
        ## removed the big title and subtitle  here
        ##title = TeX(title),
        theme = theme(plot.title =
                        element_text(size = 20,
                                     family = "CMU Serif",
                                     hjust = 0.5))
      ) &
      theme(legend.position = "bottom") &
      scale_y_continuous(breaks = c(5, 10, 20, 40, 60, 80, 90, 100))


    ggsave(paste0(parameter,
                  " (",
                  outcome,
                  ").pdf"),
           p.final,
           path = dir_plots,
           width = 21,
           height = 21.5*30/37 - 1,
           units = "cm")

  }
}

## @knitr results_plots_measures_pairs

for (parameter in c("Correlation", "Covariance")) {
  for (upper_lower in list(c("Bias", "RMSE"),
                           c("Bias", "Variance"),
                           c("Bias", "MAE"),
                           c("Positive error", "Negative error"),
                           c("Positive error", "Total error"),
                           c("Negative error", "Total error"))) {
    upper <- upper_lower[1]
    lower <- upper_lower[2]

    p.q <- plot_quadrants(d_important,
                           parameter,
                           upper,
                           lower)
      ## Saving the plot
      ggsave(
        paste0(parameter,
               " (",
               upper,
               " and ",
               lower,
               ").pdf"),
        p.q,
        path = dir_plots,
        width = 21,
        height = 30,
        units = "cm"
      )
  }
}

## @knitr results_plots_caterpillar

# Making lists of level-2 distributions, model names, and analyses --------

list_l2.dist <- c("Gaussian",
                  "Chi2")

list_model <- c("NAR",
                "Chi2AR",
                "BinAR",
                "PoDAR")

list_model.type <- c("resid.fixed",
                     "resid.random")

list_parameter <- c("Correlation",
                    "Covariance")

# Plots with 9 panels -----------------------------------------------------

for (single_parameter in list_parameter) {

cross3(list_model,
       list_l2.dist,
       list_model.type) %>%
  plyr::l_ply(function(x){
    single_model <- x[1]
    single_model_name <- ifelse(single_model == "Chiar",
                                "$\\chi^2$AR(1)",
                                paste0(single_model, "(1)"))

    single_l2.dist <- x[2]
    single_l2.dist_name <- ifelse(single_l2.dist == "Chi2",
                                  "$\\chi^2$",
                                  single_l2.dist)

    model_type <- x[3]

    caterpillar_single <- plot_caterpillar(d_abridged,
                                            l2.dist_ = single_l2.dist,
                                            Model_ = single_model,
                                            parameter_ = single_parameter,
                                            analysis.type = model_type)

    p.final <- caterpillar_single +
      plot_layout(guides = "collect") +
      plot_annotation(
        theme = theme(plot.title =
                        element_text(size = 12,
                                     family = "CMU Serif",
                                     hjust = 0.5),
                      plot.subtitle =
                        element_text(size = 12,
                                     family = "CMU Serif",
                                     hjust = 0.5))
      ) &
      theme(legend.position = "bottom")


    ggsave(paste0("caterpillar-",
                  single_parameter,
                  "-",
                  single_model,
                  "-",
                  single_l2.dist,
                  "-",
                  model_type,
                  ".pdf"),
           p.final,
           path = dir_plots,
           width = 18,
           height = 11,
           units = "cm")

  })
}

# Plots with 18 panels ----------------------------------------------------


cross3(list_model,
       list_l2.dist,
       list_parameter) %>%
  plyr::l_ply(function(x){

    single_model <- x[1]
    single_model_name <- ifelse(single_model == "Chiar",
                                "$\\chi^2$AR(1)",
                                paste0(single_model, "(1)"))

    single_l2.dist <- x[2]
    single_l2.dist_name <- ifelse(single_l2.dist == "Chi2",
                                  "$\\chi^2$",
                                  single_l2.dist)

    single_parameter <- x[3]
    caterpillar_fixed <- plot_caterpillar(d_abridged,
                                          l2.dist_ = single_l2.dist,
                                          Model_ = single_model,
                                          parameter_ = single_parameter,
                                          analysis.type = "resid.fixed")
    caterpillar_random <- plot_caterpillar(d_abridged,
                                           l2.dist_ = single_l2.dist,
                                           Model_ = single_model,
                                           parameter_ = single_parameter,
                                           analysis.type = "resid.random")

    p.patchwork <- caterpillar_fixed / caterpillar_random

    p.final <- p.patchwork +
      plot_layout(guides = "collect") +
      plot_annotation(
        theme = theme(plot.title =
                        element_text(size = 12,
                                     family = "CMU Serif",
                                     hjust = 0.5),
                      plot.subtitle =
                        element_text(size = 12,
                                     family = "CMU Serif",
                                     hjust = 0.5))
      ) &
      theme(legend.position = "bottom")


    ggsave(paste0("caterpillar-",
                  single_parameter,
                  "-",
                  single_model,
                  "-",
                  single_l2.dist,
                  "-",
                  "both",
                  ".pdf"),
           p.final,
           path = dir_plots,
           width = 18,
           height = 18,
           units = "cm")

  })


## @knitr make_dgm_df

make_dgm_df <- function(obj.1,
                        obj.2,
                        obj.3 = NULL){

  d.1 <- data.frame(t = 1:length(obj.1$x),
                    value = obj.1$x,
                    obj.id = obj.1$Model.Description.Short)

  d.2 <- data.frame(t = 1:length(obj.2$x),
                    value = obj.2$x,
                    obj.id = obj.2$Model.Description.Short)

  d <- rbind(d.1, d.2)

  if(!is.null(obj.3)){
    d.3 <- data.frame(t = 1:length(obj.3$x),
                      value = obj.3$x,
                      obj.id = obj.3$Model.Description.Short)
    d <- rbind(d, d.3)
  }

  d <- d %>%
    group_by(obj.id) %>%
    mutate(Mean = mean(value))

  return(d)
}

## @knitr additional_plots_dgms

ts.length <- 5000


d.3.nar <- make_dgm_df(
  dgm_nar(
    phi = 0.4,
    var.resid = 20,
    T = ts.length,
    Mean = 85,
    seed = 1 + 9
  ),
  dgm_nar(
    phi = 0.8,
    var.resid = 20,
    T = ts.length,
    Mean = 55,
    seed = 2 + 4
  ),
  dgm_nar(
    phi = 0.4,
    var.resid = 47,
    T = ts.length,
    Mean = 20,
    seed = 3 + 2
  )
)


profile_nar <- plot_dgm.profile(d.3.nar,
                                brewer.pal(name = "YlOrBr", n = 9)[c(5, 7, 9)],
                                "$AR(1)$")


d.3.chiar <- make_dgm_df(
  dgm_generator(
    Model = "ChiAR",
    phi = 0.4,
    nu = 25,
    T = ts.length,
    seed = 1 + 1
  ),
  dgm_generator(
    Model = "ChiAR",
    phi = 0.7,
    nu = 5,
    T = ts.length,
    seed = 2 + 5
  ),
  dgm_generator(
    Model = "ChiAR",
    phi = 0.4,
    nu = 1,
    T = ts.length,
    seed = 3 + 1
  )
)


profile_chiar <- plot_dgm.profile(d.3.chiar,
                                  brewer.pal(name = "GnBu", n = 9)[c(5, 7, 9)],
                                  "$\\chi^2AR(1)$")


d.3.binar <- make_dgm_df(
  dgm_generator(
    Model = "BinAR",
    k = 7,
    alpha  = 0.85,
    phi = 0.45,
    T = ts.length,
    seed = 1 + 5
  ),
  dgm_generator(
    Model = "BinAR",
    k = 7,
    alpha = 0.85,
    phi = 0.7,
    T = ts.length,
    seed = 2 + 6
  ),
  dgm_generator(
    Model = "BinAR",
    k = 7,
    alpha = 0.5,
    phi = 0.45,
    T = ts.length,
    seed = 3
  )
)

profile_binar <- plot_dgm.profile(d.3.binar,
                                  brewer.pal(name = "YlGn", n = 9)[c(5, 7, 9)],
                                  "$BinAR(1)$")


d.3.podar <- make_dgm_df(
  dgm_generator(
    Model = "PoDAR",
    tau = 0.4,
    Mean = 40,
    T = ts.length,
    seed = 1 + 4
  ),
  dgm_generator(
    Model = "PoDAR",
    tau = 0.8,
    Mean = 10,
    T = ts.length,
    seed = 2 + 4
  ),
  dgm_generator(
    Model = "PoDAR",
    tau = 0.4,
    Mean = 1.,
    T = ts.length,
    seed = 3 + 2
  )
)


profile_podar <- plot_dgm.profile(d.3.podar,
                                  brewer.pal(name = "BuPu", n = 9)[c(5, 7, 9)],
                                  "$\\PoDAR(1)$")

## @knitr additional_plots_dgms_save


ggsave(
  "Profile NAR.pdf",
  profile_nar,
  path = dir_plots,
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)

ggsave(
  "Profile Chi2AR.pdf",
  profile_chiar,
  path = dir_plots,
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)


ggsave(
  "Profile BinAR.pdf",
  profile_binar,
  path = dir_plots,
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)


ggsave(
  "Profile PoDAR.pdf",
  profile_podar,
  path = dir_plots,
  height = 10 * 1.8,
  width = 21 * 2,
  units = "cm"
)


ggsave(
  "Profile four DGMs.pdf",
  profile_nar /
    profile_chiar /
    profile_binar /
    profile_podar,
  path = dir_plots,
  height = 10 * 1.8 * 4,
  width = 21 * 2,
  units = "cm"
)


ggsave(
  "Profile alternative DGMs.pdf",
  profile_chiar /
    profile_binar /
    profile_podar,
  path = dir_plots,
  height = 10 * 1.8 * 3,
  width = 21 * 2,
  units = "cm"
)


## @knitr additional_plots_dataset_profiles

# Reading simulated datasets and saving in lists --------------------------

dgm_datasets_gaussian.means <- list(
  nar = "sim_uSeed-13330_l2.dist-Gaussian_Model-NAR",
  chiar = "sim_uSeed-13297_l2.dist-Gaussian_Model-Chi2AR",
  binar = "sim_uSeed-13286_l2.dist-Gaussian_Model-BinAR",
  podar = "sim_uSeed-13319_l2.dist-Gaussian_Model-PoDAR"
) %>%
  lapply(function(x) read_rds(paste0("simulation-files/fit-files/",
                                     x,
                                     "_N-100_phi-0.4_T-100_sim.Seed-0_Rep-1.rds")))


dgm_datasets_chi2.means <- list(
  nar = "sim_uSeed-13337_l2.dist-Chi2_Model-NAR",
  chiar = "sim_uSeed-13304_l2.dist-Chi2_Model-Chi2AR",
  binar = "sim_uSeed-13293_l2.dist-Chi2_Model-BinAR",
  podar = "sim_uSeed-13326_l2.dist-Chi2_Model-PoDAR"
) %>%
  lapply(function(x) read_rds(paste0("simulation-files/fit-files/",
                                     x,
                                     "_N-100_phi-0.4_T-100_sim.Seed-0_Rep-1.rds")))



# Plot profiles per DGM ---------------------------------------------------

purrr::pwalk(
  list(
    dgm_datasets_gaussian.means,
    dgm_datasets_chi2.means,
    names(dgm_datasets_gaussian.means)
  ),
  ~ save_dataset_profile(..1,
                         ..2,
                         paste0("profiles-dataset-", ..3))

)

