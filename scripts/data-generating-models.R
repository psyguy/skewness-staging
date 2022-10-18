# Data-generating models specifications -----------------------------------

## @knitr dgm_nar

dgm_nar <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  ## setting default seed if not given
  if(is.null(pa$phi)) pa$phi <- 0.4
  if(is.null(pa$Mean)) pa$Mean <- 50
  if(is.null(pa$var.marginal)) pa$var.marginal <- 4
  if(is.null(pa$var.resid)) pa$var.resid <- pa$var.marginal * (1 - pa$phi ^ 2)
  if(is.null(pa$k)) pa$k <- 100
  if(is.null(pa$T)) pa$T <- 100
  if(is.null(pa$seed)) pa$seed <- 0

  ## Making sure var.marginal and var.resid are correctly related
  pa$var.marginal <- pa$var.resid / (1 - pa$phi ^ 2)

  ### first make time series centered around zero
  ## drawing the first sample x_cent_1
  set.seed(pa$seed)
  x_cent <- rep(NA, pa$T)
  x_cent[1] <- rnorm(n = 1,
                     mean = 0,
                     sd = sqrt(pa$var.marginal)
  )

  ## making the rest of the centered time series
  for (t in 2:pa$T){

    x_cent[t] <- pa$phi*x_cent[t-1] + rnorm(n = 1,
                                            mean = 0,
                                            sd = sqrt(pa$var.resid)
    )
  }
  ## adding the mean to the centered time series
  x <- x_cent + pa$Mean

  ## quick output of raw time series without book-keeping variables/parameters
  if(!is.null(pa$only.ts))
    if(pa$only.ts==TRUE) return(x)

  Empirical.Parameters = list(Mean = mean(x),
                              Variance = var(x),
                              Skewness = moments::skewness(x),
                              AR = acf(x, lag.max = 1, plot = FALSE)$acf[2]
  )

  ## making a LaTeX-ready list description of the model
  Model.Description <- paste0("\\mu = ",
                              round(pa$Mean,2),
                              "(",
                              round(Empirical.Parameters$Mean,2),
                              ")",
                              ",\\; \\gamma = ",
                              0,
                              "(",
                              round(Empirical.Parameters$Skewness,2),
                              ")",
                              ",\\; \\phi = ",
                              round(pa$phi,2),
                              "(",
                              round(Empirical.Parameters$AR,2),
                              ")",
                              ",\\; \\sigma^2_{marginal} = ",
                              round(pa$var.marginal,2),
                              "(",
                              round(Empirical.Parameters$Variance,2),
                              ")",
                              ",\\; T = ",
                              pa$T,
                              "$")

  Model.Description.Short <- paste0("$\\mu = ",
                                    round(pa$Mean, 2),
                                    ",\\; \\sigma^2_{\\epsilon} = ",
                                    round(pa$var.resid, 2),
                                    ",\\; \\phi = ",
                                    round(pa$phi, 2),
                                    "\\; \\rightarrow",
                                    "\\; \\sigma^2 = ",
                                    round(Empirical.Parameters$Variance,2),
                                    ",\\; \\gamma = ",
                                    round(Empirical.Parameters$Skewness,2),
                                    "$")


  ## making the output object
  output <- list(x = x,
                 Model.Description = Model.Description,
                 Model.Description.Short = Model.Description.Short,
                 Model.Parameters = pa,
                 Empirical.Parameters = Empirical.Parameters)

  return(output)
}

## @knitr dgm_chiar

dgm_chiar <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  ## setting default seed if not given
  if(is.null(pa$phi)) pa$phi <- 0.4
  if(is.null(pa$nu)) pa$nu <- 3
  if(is.null(pa$c)) pa$c <- 0
  if(is.null(pa$k)) pa$k <- 100
  if(is.null(pa$T)) pa$T <- 100
  if(is.null(pa$seed)) pa$seed <- 0

  ## drawing the first sample x_1
  set.seed(pa$seed)
  x <- rep(NA, pa$T)
  x[1] <- rchisq(n = 1,
                 df = pa$Mean)

  ## making the rest of the time series
  for (t in 2:pa$T){

    x[t] <- pa$c + pa$phi*x[t-1] + rchisq(n = 1,
                                          df = pa$nu)
  }

  ## quick output of raw time series without book-keeping variables/parameters
  if(!is.null(pa$only.ts))
    if(pa$only.ts==TRUE) return(x)

  Empirical.Parameters = list(Mean = mean(x),
                              Variance = var(x),
                              Skewness = moments::skewness(x),
                              AR = acf(x, lag.max = 1, plot = FALSE)$acf[2]
  )

  ## making a LaTeX-ready list description of the model
  Model.Description <- paste0("$\\chi^2AR(1):", # \\; with",
                              "\\; \\mu = ",
                              round(pa$Mean,2),
                              "(",
                              round(Empirical.Parameters$Mean,2),
                              ")",
                              ",\\; \\gamma = ",
                              round(pa$Skewness,2),
                              "(",
                              round(Empirical.Parameters$Skewness,2),
                              ")",
                              ",\\; \\phi = ",
                              round(pa$phi,2),
                              "(",
                              round(Empirical.Parameters$AR,2),
                              ")",
                              ",\\; \\nu = ",
                              round(pa$nu,2),
                              ",\\; c = ",
                              round(pa$c,3),
                              ",\\; T = ",
                              pa$T,
                              "$")

  Model.Description.Short <- paste0("$c = ",
                                    round(pa$c,2),
                                    ",\\; \\nu = ",
                                    round(pa$nu, 2),
                                    ",\\; \\phi = ",
                                    round(pa$phi, 2),
                                    "\\; \\rightarrow",
                                    "\\; \\mu = ",
                                    round(Empirical.Parameters$Mean, 2),
                                    ",\\; \\sigma^2 = ",
                                    round(Empirical.Parameters$Variance,2),
                                    ",\\; \\gamma = ",
                                    round(Empirical.Parameters$Skewness,2),
                                    "$")
  ## making the output object
  output <- list(x = x,
                 Model.Description = Model.Description,
                 Model.Description.Short = Model.Description.Short,
                 Model.Parameters = pa,
                 Empirical.Parameters = Empirical.Parameters)

  return(output)
}

## @knitr dgm_binar

dgm_binar <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  ## setting default seed if not given
  if(is.null(pa$alpha)) pa$alpha <- 0.5
  if(is.null(pa$beta)) pa$beta <- 0.4
  if(is.null(pa$k)) pa$k <- 10
  if(is.null(pa$T)) pa$T <- 100
  if(is.null(pa$seed)) pa$seed <- 0

  ## making other parameters
  pa$rho <- pa$alpha - pa$beta
  pa$theta <- pa$beta/(1-pa$rho)


  ## drawing the first sample x_1
  set.seed(pa$seed)
  x <- rep(NA, pa$T)
  x[1] <- rbinom(n = 1,
                 size = pa$k ,
                 prob = pa$theta)

  ## making the rest of the time series
  for (t in 2:pa$T){
    S_t <- rbinom(n = 1,
                  size = x[t-1],
                  prob = pa$alpha)
    R_t <- rbinom(n = 1,
                  size = pa$k - x[t-1],
                  prob = pa$beta)
    x[t] <- S_t + R_t
  }

  ## quick output of raw time series without book-keeping variables/parameters
  if(!is.null(pa$only.ts))
    if(pa$only.ts==TRUE) return(x)

  Empirical.Parameters = list(Mean = mean(x),
                              Variance = var(x),
                              Skewness = moments::skewness(x),
                              AR = acf(x, lag.max = 1, plot = FALSE)$acf[2]
  )

  ## making a LaTeX-ready list description of the model
  Model.Description <- paste0("$BinAR(1):",
                              "\\; \\mu = ",
                              round(pa$Mean,2),
                              "(",
                              round(Empirical.Parameters$Mean,2),
                              ")",
                              ",\\; \\gamma = ",
                              round(pa$Skewness,2),
                              "(",
                              round(Empirical.Parameters$Skewness,2),
                              ")",
                              ",\\; \\rho = ",
                              round(pa$rho,3),
                              "(",
                              round(Empirical.Parameters$AR,2),
                              ")",
                              ",\\; \\alpha = ",
                              round(pa$alpha,2),
                              ",\\; \\beta = ",
                              round(pa$beta,2),
                              ",\\; \\theta = ",
                              round(pa$theta,3),
                              ",\\; T = ",
                              pa$T,
                              "$")

  Model.Description.Short <- paste0("$k =",
                                    round(pa$k,2),
                                    ",\\; \\alpha = ",
                                    round(pa$alpha,2),
                                    ",\\; \\beta = ",
                                    round(pa$beta,2),
                                    "\\; \\rightarrow",
                                    "\\; \\rho = ",
                                    round(Empirical.Parameters$AR, 2),
                                    ",\\; \\mu = ",
                                    round(Empirical.Parameters$Mean, 2),
                                    ",\\; \\sigma^2 = ",
                                    round(Empirical.Parameters$Variance,2),
                                    ",\\; \\gamma = ",
                                    round(Empirical.Parameters$Skewness,2),
                                    "$")
  ## making the output object
  output <- list(x = x,
                 Model.Description = Model.Description,
                 Model.Description.Short = Model.Description.Short,
                 Model.Parameters = pa,
                 Empirical.Parameters = Empirical.Parameters)

  return(output)
}

## @knitr dgm_podar

dgm_podar <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  ## setting default seed if not given
  if(is.null(pa$tau)) pa$tau <- 0.7
  if(is.null(pa$lambda)) pa$lambda <- 0.5
  if(is.null(pa$k)) pa$k <- 6
  if(is.null(pa$T)) pa$T <- 100
  if(is.null(pa$seed)) pa$seed <- 0


  ## drawing the first sample x_1
  set.seed(pa$seed)
  x <- rep(NA, pa$T)
  x[1] <- rpois(n = 1,
                lambda = pa$lambda)

  ## making the rest of the time series
  for (t in 2:pa$T){
    V_t <- rbinom(n = 1,
                  size = 1,
                  prob = pa$tau)
    Z_t <- rpois(n = 1,
                 lambda = pa$lambda)
    x[t] <- V_t*x[t-1] + (1-V_t)*Z_t
  }

  ## quick output of raw time series without book-keeping variables/parameters
  if(!is.null(pa$only.ts))
    if(pa$only.ts==TRUE) return(x)

  Empirical.Parameters = list(Mean = mean(x),
                              Variance = var(x),
                              Skewness = moments::skewness(x),
                              AR = acf(x, lag.max = 1, plot = FALSE)$acf[2]
  )

  ## making a LaTeX-ready list description of the model
  Model.Description <- paste0("$PoDAR(1):",
                              "\\; \\mu = ",
                              round(pa$Mean,2),
                              "(",
                              round(Empirical.Parameters$Mean,2),
                              ")",
                              ",\\; \\gamma = ",
                              round(pa$Skewness,2),
                              "(",
                              round(Empirical.Parameters$Skewness,2),
                              ")",
                              "\\; \\tau = ",
                              round(pa$tau,3),
                              "(",
                              round(Empirical.Parameters$AR,2),
                              ")",
                              ",\\; \\lambda = ",
                              round(pa$lambda,3),
                              ",\\; T = ",
                              pa$T,
                              "$")


  Model.Description.Short <- paste0("$\\lambda = ",
                                    round(pa$lambda,3),
                                    ",\\; \\tau = ",
                                    round(pa$tau, 2),
                                    "\\; \\rightarrow",
                                    "\\; \\rho = ",
                                    round(Empirical.Parameters$AR, 2),
                                    ",\\; \\mu = ",
                                    round(Empirical.Parameters$Mean, 2),
                                    ",\\; \\sigma^2 = ",
                                    round(Empirical.Parameters$Variance,2),
                                    ",\\; \\gamma = ",
                                    round(Empirical.Parameters$Skewness,2),
                                    "$")
  ## making the output object
  output <- list(x = x,
                 Model.Description = Model.Description,
                 Model.Description.Short = Model.Description.Short,
                 Model.Parameters = pa,
                 Empirical.Parameters = Empirical.Parameters)


  return(output)
}

# General DGM wrappers ----------------------------------------------------

## @knitr dgm_parameterizer

dgm_parameterizer <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  if(is.null(pa$Model)) pa$Model <- "ChiAR(1)"
  if(is.null(pa$phi)) pa$phi <- 0.2

  ## %%%%%%%%%%%%
  ## for NAR(1)
  ## %%%%%%%%%%%%

  if(tolower(pa$Model) == "nar(1)" | tolower(pa$Model) == "nar"){

    if(is.null(pa$k)) pa$k <- 100

    ## Calculating model parameters
    ## if mean is given
    if (!is.null(pa$Mean)) {
      # then from the mean formula
      pa$c <- pa$Mean * (1 - pa$phi)
    }
    ## if mean is not given, we get it from the intercept
    else{
      # set intercept to zero, if already not defined
      if(is.null(pa$c)) pa$c <- 0
      # then from the mean formula
      pa$Mean <- pa$c / (1 - pa$phi)
    }

    ## Now we certainly have the mean calculated

    ## The within-person (marginal) variance is more interpretable. So
    ## we always calculate the residual variance based on the marginal
    ## variance, unless var.resid is known BUT var.marginal is NOT.
    if (!is.null(pa$var.resid) & is.null(pa$var.marginal)) {
      # from the marginal variance formula
      pa$var.marginal <- pa$var.resid / (1 - pa$phi ^ 2)
    }
    ## If none of them are known we first give var.marginal a default
    ## value and calculate var.resid based on that.
    if (is.null(pa$var.resid) & is.null(pa$var.marginal)) {
      # we set a default value for the marginal variance
      pa$var.marginal <- 4
    }
    ## Then, we calculate var.resid (again) with var.marginal to make
    ## sure that the marginal variance had been given more importance
    pa$var.resid <- pa$var.marginal * (1 - pa$phi ^ 2)

    ## Returning the parameter list
    return(pa)

  }


  ## %%%%%%%%%%%%
  ## for ChiAR(1)
  ## %%%%%%%%%%%%

  if(tolower(pa$Model) == "chiar(1)" | tolower(pa$Model) == "chiar" |
     tolower(pa$Model) == "chi2ar(1)" | tolower(pa$Model) == "chi2ar"){

    if(is.null(pa$k)) pa$k <- 100

    ## Calculating model parameters
    ## if mean and skewness are given
    if (!is.null(pa$Mean) & !is.null(pa$Skewness)) {
      # from the skewness formula
      pa$nu <-
        8 * (1 - pa$phi ^ 2) ^ 3 / ((pa$Skewness ^ 2) * (1 - pa$phi ^
                                                           3) ^ 2)
      # then from the mean formula
      pa$c <- pa$Mean * (1 - pa$phi) - pa$nu
      # returning the parameter list
      return(pa)
    }
    ## if mean and skewness are not given at the same time, we use c
    ## and one other parameter
    else{

      # set intercept to zero, if already not defined
      if(is.null(pa$c)) pa$c <- 0

      ## if mean is given
      if (!is.null(pa$Mean)) {
        # then from the mean formula
        pa$nu <- pa$Mean * (1 - pa$phi) - pa$c
        # from the skewness formula
        pa$Skewness <-
          2 * (1 - pa$phi ^ 2) ^ 1.5 / (sqrt(pa$nu / 2) * (1 - pa$phi ^
                                                             3))
        # returning the parameter list
        return(pa)
      }

      ## if skewness is given
      if (!is.null(pa$Skewness)) {
        # from the skewness formula
        pa$nu <-
          8 * (1 - pa$phi ^ 2) ^ 3 / ((pa$Skewness ^ 2) * (1 - pa$phi ^
                                                             3) ^ 2)
        # then from the mean formula
        pa$Mean <- (pa$c + pa$nu) / (1 - pa$phi)
        # returning the parameter list
        return(pa)
      }

      ## if nu is given
      if(!is.null(pa$nu)){
        # from the mean formula
        pa$Mean <- (pa$c + pa$nu)/(1-pa$phi)
        # from the skewness formula
        pa$Skewness <-
          2 * (1 - pa$phi ^ 2) ^ 1.5 / (sqrt(pa$nu / 2) * (1 - pa$phi ^
                                                             3))
        # returning the parameter list
        return(pa)
      }


    }


  }


  ## %%%%%%%%%%%%
  ## for BinAR(1)
  ## %%%%%%%%%%%%

  if (tolower(pa$Model) == "binar(1)" | tolower(pa$Model) == "binar") {

    if (is.null(pa$k))
      pa$k <- 10

    ## Calculating model parameters

    ## if mean is given
    if (!is.null(pa$Mean)) {
      # from skewness formula: m = k*theta
      pa$theta <- pa$Mean / pa$k
      ## we then calculate skewness based on theta
      pa$Skewness <- (1 - 2 * pa$theta) / sqrt(pa$k * pa$theta * (1 - pa$theta))
      ## we then calculate beta based on theta and phi
      pa$beta <- pa$theta * (1 - pa$phi)
      # then we calculate alpha
      pa$alpha <- pa$phi + pa$beta
      # finally we calculate k*beta, which is equivalent to c
      pa$c <- pa$k * pa$beta
      # returning the parameter list
      return(pa)
    }

    ## if skewness is given
    if (!is.null(pa$Skewness)) {
      # from skewness formula: skewness = (1-2*theta)/sqrt(k*theta*(1-theta))
      # it is easier to write with ks = k*(skewness^2)
      ks <- pa$k*pa$Skewness^2
      pa$theta <- (ks + 4 -
                     sqrt(ks * (ks + 4))) / (2 * ks + 8)
      # The above formula only gives theta < 0.5, thus for negative skewness
      # we must use 1-theta instead
      if(pa$Skewness < 0) pa$theta <- 1 - pa$theta
      # we then calculate mean based on theta
      pa$Mean <- pa$k * pa$theta
      # we then calculate beta based on theta and phi
      pa$beta <- pa$theta * (1 - pa$phi)
      # then we calculate alpha
      pa$alpha <- pa$phi + pa$beta
      # finally we calculate k*beta, which is equivalent to c
      pa$c <- pa$k * pa$beta
      # returning the parameter list
      return(pa)
    }

    ## if theta is given
    if (!is.null(pa$theta)) {
      # from theta formula
      pa$beta <- pa$theta * (1 - pa$phi)
      # from beta formula
      pa$alpha <- pa$phi + pa$beta
      # we then calculate mean based on theta
      pa$Mean <- pa$k * pa$theta
      # we then calculate skewness based on theta
      pa$Skewness <- (1 - 2 * pa$theta) / sqrt(pa$k * pa$theta * (1 - pa$theta))
      # finally we calculate k*beta, which is equivalent to c
      pa$c <- pa$k * pa$beta
      # returning the parameter list
      return(pa)
    }

    ## if alpha is given
    if (!is.null(pa$alpha)) {
      # from beta formula
      pa$beta <- pa$alpha - pa$phi
      # from theta formula
      pa$theta <- pa$beta / (1 - pa$phi)
      # we then calculate mean based on theta
      pa$Mean <- pa$k * pa$theta
      # we then calculate skewness based on theta
      pa$Skewness <- (1 - 2 * pa$theta) / sqrt(pa$k * pa$theta * (1 - pa$theta))
      # finally we calculate k*beta, which is equivalent to c
      pa$c <- pa$k * pa$beta
      # returning the parameter list
      return(pa)
    }

    ## if beta is given
    if (!is.null(pa$beta)) {
      # from beta formula
      pa$alpha <- pa$phi + pa$beta
      # from theta formula
      pa$theta <- pa$beta / (1 - pa$phi)
      # we then calculate mean based on theta
      pa$Mean <- pa$k * pa$theta
      # we then calculate skewness based on theta
      pa$Skewness <- (1 - 2 * pa$theta) / sqrt(pa$k * pa$theta * (1 - pa$theta))
      # finally we calculate k*beta, which is equivalent to c
      pa$c <- pa$k * pa$beta
      # returning the parameter list
      return(pa)
    }

    ## if intercept given
    if (!is.null(pa$c)) {
      # from intercept formula c = k*beta
      pa$beta <- pa$c / pa$k
      # from beta formula
      pa$alpha <- pa$phi + pa$beta
      # from theta formula
      pa$theta <- pa$beta / (1 - pa$phi)
      # we then calculate mean based on theta
      pa$Mean <- pa$k * pa$theta
      # we then calculate skewness based on theta
      pa$Skewness <- (1 - 2 * pa$theta) / sqrt(pa$k * pa$theta * (1 - pa$theta))
      # returning the parameter list
      return(pa)
    }

  }

  ## %%%%%%%%%%%%
  ## for PoDAR(1)
  ## %%%%%%%%%%%%

  if (tolower(pa$Model) == "podar(1)" | tolower(pa$Model) == "podar") {

    if (is.null(pa$k))
      pa$k <- 50

    # phi and tau are the same, then if tau is defined, it overrules phi
    if(!is.null(pa$tau)) pa$phi <- pa$tau
    # and if tau is not defined, then tau will get the value of phi
    if(is.null(pa$tau)) pa$tau <- pa$phi
    # DAR(1) intercept is zero
    pa$c <- 0

    ## Calculating model parameters

    ## if mean is given
    if (!is.null(pa$Mean)) {
      # mean of Poisson is lambda
      pa$lambda <- pa$Mean
      ## we then calculate skewness based on lambda
      pa$Skewness <- pa$lambda^(-0.5)
      # returning the parameter list
      return(pa)
    }

    ## if skewness is given
    if (!is.null(pa$Skewness)) {
      # from skewness formula: skewness = 1/sqrt(lambda)
      pa$lambda <- pa$Skewness^(-2)
      # we then calculate mean based on theta
      pa$Mean <- pa$lambda
      # returning the parameter list
      return(pa)
    }

    ## if lambda is given
    if (!is.null(pa$lambda)) {
      # mean of Poisson is lambda
      pa$Mean <- pa$lambda
      ## we then calculate skewness based on lambda
      pa$Skewness <- pa$lambda^(-0.5)
      # returning the parameter list
      return(pa)
    }
  }

}

## @knitr dgm_generator

dgm_generator <- function(...){

  pa <- list(...)

  if(is.list(pa$pa)) pa <- pa$pa

  ## setting default seed if not given
  if(is.null(pa$Model)) pa$Model <- "ChiAR(1)"
  if(is.null(pa$phi)) pa$phi <- 0.2
  # if no model parameter is given, then mean is set to a default 5
  if (is.null(pa$Mean) &
      is.null(pa$Skewness) &
      is.null(pa$c) &
      is.null(pa$nu) &
      is.null(pa$alpha) &
      is.null(pa$beta) &
      is.null(pa$theta)
  ) pa$Mean <- 5
  # if(is.null(pa$Mean)) pa$Variance <- 10
  # if(is.null(pa$Skewness)) pa$Skewness <- 3
  if(is.null(pa$T)) pa$T <- 100
  if(is.null(pa$seed)) pa$seed <- 0

  if(is.null(pa$only.ts)) pa$only.ts <- FALSE

  ## If you set `pa$only.ts` parameter as TRUE, the dgm_ functions produce only
  ## the raw time series (a single numeric vector) which is way lighter and way
  ## faster:

  # pa$only.ts <- TRUE


  ## calculating model parameters
  pa <- dgm_parameterizer(pa = pa)

  ### making models

  ## %%%%%%%%%%%%
  ## NAR(1)
  ## %%%%%%%%%%%%
  if(tolower(pa$Model) == "nar(1)" | tolower(pa$Model) == "nar"){
    # default maximum scale value
    if(is.null(pa$k)) pa$k <- 100
    ## %% Generating the data
    o <- dgm_nar(pa = pa)
  }

  ## %%%%%%%%%%%%
  ## ChiAR(1)
  ## %%%%%%%%%%%%
  if(tolower(pa$Model) == "chiar(1)" | tolower(pa$Model) == "chiar" |
     tolower(pa$Model) == "chi2ar(1)" | tolower(pa$Model) == "chi2ar"){
    # default maximum scale value
    if(is.null(pa$k)) pa$k <- 100
    ## %% Generating the data
    o <- dgm_chiar(pa = pa)
  }

  ## %%%%%%%%%%%%
  ## BinAR(1)
  ## %%%%%%%%%%%%
  if(tolower(pa$Model) == "binar(1)" | tolower(pa$Model) == "binar"){
    # default maximum scale value
    if(is.null(pa$k)) pa$k <- 10
    ## %% Generating the data
    o <- dgm_binar(pa = pa)
  }

  ## %%%%%%%%%%%%
  ## DAR(1)
  ## %%%%%%%%%%%%
  if(tolower(pa$Model) == "dar(1)" | tolower(pa$Model) == "dar"){
    # default maximum scale value
    if(is.null(pa$k)) pa$k <- 10
    ## %% Generating the data
    o <- dgm_dar(pa = pa)
  }

  ## %%%%%%%%%%%%
  ## PoDAR(1)
  ## %%%%%%%%%%%%
  if(tolower(pa$Model) == "podar(1)" | tolower(pa$Model) == "podar"){
    # default maximum scale value
    if(is.null(pa$k)) pa$k <- 100
    ## %% Generating the data
    o <- dgm_podar(pa = pa)
  }

  ## Also allow a data frame output
  # if(!is.null(p$as.dataframe) & p$as.dataframe){
  #
  # }

  return(o)
}


# Dataset generation ------------------------------------------------------

## @knitr dgm_make.sample

dgm_make.sample <- function(Model = "ChiAR(1)",
                            Means = rnorm(100, 5, 3),
                            T = 100,
                            phi = 0.4,
                            seeds = NULL
){

  N <- length(Means)
  df <- data.frame(subject = rep(1:N, each = T),
                   t = rep(1:T, times = N),
                   x = rep(NA, N*T))
  if(is.null(seeds)) seeds <- 1000*N*Means
  if(length(seeds)<=1) seeds <- seeds + 1000*N*Means

  for(s in 1:N){
    x <- dgm_generator(
      Model = Model,
      only.ts = TRUE,
      T = T,
      phi = phi,
      Mean = Means[s],
      # Skewness = Skewness,
      seed = seeds[s])

    df$x[((s-1)*T + 1):(s*T)] <- x

  }

  return(df)
}

## @knitr make_datasets

make_datasets <- function(Model = "DAR",
                          T = 100,
                          N = 100,
                          phi = 0.4,
                          l2.distribution = "Gaussian",
                          seed = 0) {
  # save global seed of the global env and set it back before leaving
  seed.old <- .Random.seed
  on.exit({
    .Random.seed <<- seed.old
  })
  set.seed(seed)

  if (tolower(Model) == "chiar" | tolower(Model) == "chi2ar") {
    model.name <- "Chi2AR"
    lev2.Mean <- 10
    lev2.Variance <- 10
    chi2.df <- 5
  }
  if (tolower(Model) == "binar") {
    model.name <- "BinAR"
    lev2.Mean <- 2
    lev2.Variance <- 1
    chi2.df <- 2.9
  }
  if (tolower(Model) == "dar") {
    model.name <- "DAR"
    lev2.Mean <- 2
    lev2.Variance <- 1
    chi2.df <- 2.9
  }

  if (tolower(Model) == "podar") {
    model.name <- "PoDAR"
    lev2.Mean <- 4
    lev2.Variance <- 4
    chi2.df <- 1.5
  }

  if (tolower(Model) == "nar") {
    model.name <- "NAR"
    lev2.Mean <- 50
    lev2.Variance <- 4
    chi2.df <- 2
  }

  # sampling within-person mean from level 2 distribution
  if (l2.distribution == "Gaussian")
    Means <- rnorm(2 * N, lev2.Mean, sqrt(lev2.Variance))
  if (l2.distribution == "Chi2")
    Means <- rchisq(2 * N, chi2.df)

  # removing out-of-bounds means
  Means[Means < 0] <- NA
  Means[Means > 100] <- NA
  if (model.name == "BinAR" |
      model.name == "DAR")
    Means[Means > 10] <- NA

  # keeping N samples from the in-bound means
  Means <- Means %>% na.omit() %>% sample(N)


  sample_df <- dgm_make.sample(
    Model = Model,
    Means = Means,
    T = T,
    phi = phi,
    seeds = NULL
  )
  return(sample_df)
}
