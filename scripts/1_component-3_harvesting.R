## @knitr fit_extract

fit_extract <- function(rds.file){

  m <- readRDS(rds.file)
  # Extract all parameter estimates
  est.par <- m[["fit.Dataset"]][["results"]][["parameters"]]
  # Extract .out output text (for factor scores)
  mt <- m[["fit.Dataset"]][["results"]][["output"]]

  # make empty dataframe for NAs
  empty <- data.frame(matrix(ncol = 19, nrow = 0))
  colnames(empty) <- c("uSeed", "type", "l2.dist", "Model", "N", "phi", "T",
                       "Rep", "standardization", "est", "posterior_sd", "pval",
                       "lower_2.5ci", "upper_2.5ci", "sig", "BetweenWithin",
                       "param.name", "fit.ElapsedTime", "fit.File")

  if(length(est.par) == 0) return(empty)
  if(is.null(est.par[["unstandardized"]])) return(empty)
  if(is.null(est.par[["stdyx.standardized"]])) return(empty)

  unstd <- est.par[["unstandardized"]] %>%
    mutate(param.name = paste(paramHeader,
                              param,
                              sep = ".")
    ) %>%
    select(-paramHeader:-param) %>%
    mutate(standardization = "unstd",
           .before = est)
  stdyx <- est.par[["stdyx.standardized"]] %>%
    mutate(param.name = paste(paramHeader,
                              param,
                              sep = ".")
    ) %>%
    select(-paramHeader:-param) %>%
    mutate(standardization = "stdyx",
           .before = est)

  # ## Extracting factor scores
  #
  # mt[mt == ""] <- NA
  # mt <- mt %>% na.omit()
  # mt %>% head
  #
  # phi.start <- grep("Results for Factor PHI", mt) + 2
  # phi.end <- grep("Results for X", mt) - 1
  #
  # x.start <- grep("Results for X", mt) + 2
  # x.end <- grep("TECHNICAL 1 OUTPUT", mt) - 1
  #
  # fac_phi <- read.table(text = mt[phi.start:phi.end], header = FALSE, fill = TRUE)
  # fac_x <- read.table(text = mt[x.start:x.end], header = FALSE, fill = TRUE)
  #
  # nc <- ncol(fac_phi)
  #
  # f_phi_est <- c()
  # f_phi_cluster <- c()
  # for(cc in 1:(nc/3)){
  #   f_phi_est <- c(f_phi_est, fac_phi[, 3*cc]) %>% na.omit()
  #   f_phi_cluster <- c(f_phi_cluster, fac_phi[, 3*cc - 1]) %>% na.omit()
  # }
  # fac_phi_id.sorted <- data.frame(est = f_phi_est,
  #                                 id = f_phi_cluster) %>%
  #   arrange(id) %>%
  #   pull(est)
  #
  # f_x_est <- c()
  # f_x_cluster <- c()
  # for(cc in 1:(nc/3)){
  #   f_x_est <- c(f_x_est, fac_x[, 3*cc]) %>% na.omit()
  #   f_x_cluster <- c(f_x_cluster, fac_x[, 3*cc - 1]) %>% na.omit()
  # }
  # fac_x_id.sorted <- data.frame(est = f_x_est,
  #                               id = f_x_cluster) %>%
  #   arrange(id) %>%
  #   pull(est)


  m$fit.Dataset <- NULL

  res <- unstd %>%
    rbind(stdyx) %>%
    mutate(fit.ElapsedTime = as.numeric(m[["fit.EndTime"]] - m[["fit.StartTime"]],
                                       units = "mins"),
    fit.File = gsub(".*/", "", m$fit.File)) %>%
    mutate(uSeed = m$uSeed,
           type = m$type %>% as.character(),
           l2.dist = m$l2.dist,
           Model = m$Model,
           N = m$N,
           phi = m$phi,
           T = m$T,
           Rep = m$Rep,
           # given.cor = cor(m$given.Means, m$given.Phis),
           # given.cov = cov(m$given.Means, m$given.Phis),
           # factor.cor = cor(fac_x_id.sorted, fac_phi_id.sorted),
           # factor.cov = cov(fac_x_id.sorted, fac_phi_id.sorted),
           .before = standardization
    )

  return(res)

}

