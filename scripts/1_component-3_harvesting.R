## @knitr fit_extract

fit_extract <- function(rds.file){

  m <- readRDS(rds.file)
  est.par <- m[["fit.Dataset"]][["results"]][["parameters"]]

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

  m$fit.Dataset <- NULL

  res <- unstd %>%
    rbind(stdyx) %>%
    mutate(fit.ElapsedTime = (difftime(m[["fit.EndTime"]],
                                       m[["fit.StartTime"]],
                                       units="mins")
    ) %>%
      as.numeric(),
    fit.File = gsub(".*/", "", m$fit.File)) %>%
    mutate(uSeed = m$uSeed,
           type = m$type,
           l2.dist = m$l2.dist,
           Model = m$Model,
           N = m$N,
           phi = m$phi,
           T = m$T,
           Rep = m$Rep,
           .before = standardization
    )

  return(res)

}
