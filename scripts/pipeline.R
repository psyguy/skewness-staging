## @knitr make_sim_refs

make_sim_refs <-
  function(conditions = list(T = c(30, 100),
                             N = c(100),
                             Model = c("BinAR", "Chi2AR", "DAR"),
                             l2.dist = c("Gaussian", "Chi2"),
                             phi = c(0.4)
  ),
  Reps = 100,
  seed = 0,
  save.directory = "self-sim",
  save.refs.filename = paste("sim_refs",
                             Sys.Date()),
  prefix.sim.datasets = "sim_"){

    dir.create(save.directory, showWarnings = FALSE)

    # allow Reps to be used as a vector of indexes
    if(length(Reps)<2) Reps = seq_len(Reps)


    # save global seed of the global env and set it back before leaving
    seed.old <- .Random.seed
    on.exit({
      .Random.seed <<- seed.old
    })

    # sorting conditions alphabetically
    conditions <- conditions[order(names(conditions))]

    conditions$sim.Seed <- seed
    conditions$Rep <- Reps

    n.conditions <- length(conditions)

    # making the first columns of the data frame
    d <- conditions %>%
      expand.grid(stringsAsFactors = TRUE)

    # transforming factors to numerics
    d.numeric <- d
    factor.columns <- sapply(d.numeric, is.factor)
    d.numeric[factor.columns] <-
      sapply(d.numeric[factor.columns], as.numeric)

    # getting rid of non-integers and Rep
    d.integer <- d.numeric[, -n.conditions] %>%
      apply(2,
            function(x)
              x %>%
              as.character() %>%
              gsub("\\.", "", .) %>%
              as.numeric())

    # to make unique seeds, we must sum conditions weighted by prime numbers
    # but the primes must not be among prime factors of conditions

    # conditions prime factors
    cpfs <- d.integer %>%
      unlist() %>%
      as.numeric() %>%
      unique() %>%
      primes::prime_factors() %>%
      unlist() %>%
      unique()

    primes.seq <- c(cpfs,
                    primes::generate_n_primes(ncol(d.integer) + length(cpfs)))

    primes.seq <-
      primes.seq[!(duplicated(primes.seq) |
                     duplicated(primes.seq, fromLast = TRUE))]

    uSeed <- d.integer %*% primes.seq %>%
      as.character() %>%
      paste0(as.character(d.numeric$Rep),
             .)

    d.headers <- d %>%
      mutate(
        uSeed = uSeed,
        sim.Path = save.directory,
        sim.File = NA
      ) %>%
      relocate(uSeed, .before = 1)

    for (r in 1:nrow(d.headers)) {
      only.headers <- d.headers %>%
        select(-sim.Path:-sim.File) %>%
        colnames()
      r.values <- d.headers[r, only.headers]
      factor.columns <- sapply(r.values, is.factor)
      r.values[factor.columns] <-
        sapply(r.values[factor.columns], as.character)
      d.headers[r, "sim.File"] <- only.headers %>%
        paste0("-") %>%
        paste0(r.values) %>%
        paste(collapse = "_") %>%
        paste0(".rds") %>%
        paste0(prefix.sim.datasets, .) # prefix for simulated datasets
    }

    d <- d.headers

    # getting rid of factors
    factor.columns <- sapply(d, is.factor)
    d[factor.columns] <-
      sapply(d[factor.columns], as.character)

    # making uSeed numeric
    d$uSeed <- d$uSeed %>% as.numeric()

    # Save the references data frame to a file
    if(!is.null(save.refs.filename)){
      saveRDS(d,
              file = here::here(save.directory, paste0(save.refs.filename,
                                                       ".rds"))
      )
      write.csv(d,
                file = here::here(save.directory, paste0(save.refs.filename,
                                                         ".csv")),
                row.names = FALSE)
    }

    return(d)
  }

## @knitr make_fit_refs

make_fit_refs <-
  function(sim_refs,
           hyperparameters = list(iter = c(2000),
                                  thin = c(5),
                                  type = c("resid.random",
                                           "resid.fixed")),
           save.directory = "self-sim",
           save.refs.filename = paste("fit_refs",
                                      Sys.Date()),
           prefix.fit.datasets = "fit_"){

    dir.create(save.directory, showWarnings = FALSE)


    h <- hyperparameters %>%
      expand.grid(stringsAsFactors = TRUE)

    dd <- sim_refs[rep(1:nrow(sim_refs), times = nrow(h)),]
    rownames(dd) <- NULL
    hh <- h[rep(1:nrow(h), each = nrow(sim_refs)),]
    rownames(hh) <- NULL

    d <- cbind(dd,hh)
    d <- d %>%
      mutate(fit.Path = save.directory,
             fit.File = NA
      )


    only.headers <- names(hyperparameters) %>%
      c("uSeed",
        "N",
        "T",
        .,
        "Rep")
    only.headers <- only.headers[!(only.headers %in% c("iter", "thin"))]

    for (r in 1:nrow(d)) {
      r.values <-d[r, only.headers]
      factor.columns <- sapply(r.values, is.factor)
      r.values[factor.columns] <-
        sapply(r.values[factor.columns], as.character)
      d[r, "fit.File"] <- only.headers %>%
        paste0("-") %>%
        paste0(r.values) %>%
        paste(collapse = "_") %>%
        paste0(".rds") %>%
        paste0(prefix.fit.datasets, .) # prefix for fitted datasets
      print(r)
    }

    if(!is.null(save.refs.filename)){
      saveRDS(d,
              file = here::here(save.directory,
                                paste0(save.refs.filename,
                                       ".rds"))
      )
      write.csv(d,
                file = here::here(save.directory,
                                  paste0(save.refs.filename,
                                         ".csv")),
                row.names = FALSE)
    }

    return(d)
  }


## @knitr do_sim_parallel


do_sim_parallel <-
  function(sim_refs,
           nClust = 48,
           save.directory = "self-sim",
           alternative.sim.Path = NULL,
           clusterLOG.filename = paste0("sim_clusterLOG_",
                                        Sys.Date(),
                                        ".txt"),
           sleeptime = 1 # seconds to wait before running clusters
  ){

    cl <- snow::makeSOCKcluster(nClust,
                                outfile = here::here(save.directory,
                                                     clusterLOG.filename))
    debug <- TRUE

    d <- sim_refs


    ## Start clusters:

    snow::clusterExport(cl,
                        c("d",
                          #"make_datasets",
                          "alternative.sim.Path",
                          "debug"),
                        envir = environment())
    t.snow <- snow::snow.time({

      snow::clusterApplyLB(cl = cl,
                           seq_len(nrow(d)),
                           function(i) {

                             Sys.sleep(sleeptime*(i %% nClust))

                             source(here::here("functions",
                                               "functions_data-generating-models.R"))
                             library(tidyverse)

                             # if (debug) {
                             #   cat("\nRunning iteration:",
                             #       i,
                             #       " / ",
                             #       nrow(d),
                             #       "\nTime:",
                             #       as.character(Sys.time()),
                             #       "\n")
                             #   print(d$sim.File)
                             # }


                             d_i <- as.list(d[i, ])
                             arguments <-
                               as.list(d_i[2:(length(d_i) - 4)])
                             arguments$seed <- d_i$uSeed

                             sim.Path <- ifelse(is.null(alternative.sim.Path),
                                                d_i$sim.Path,
                                                alternative.sim.Path)

                             sim.StartTime <- Sys.time()

                             # tryRes <-
                             #   try(

                             output.dataset <- do.call(make_datasets,
                                                       arguments)
                             # )
                             # if (is(tryRes,"try-error")){
                             #   if (debug){
                             #     browser()
                             #   }
                             #   return(list(error = TRUE, errorMessage = as.character(tryRes), id = d$id[i]))
                             # }
                             d_i$sim.Dataset <- output.dataset
                             d_i$sim.StartTime <- sim.StartTime
                             d_i$sim.EndTime <- Sys.time()
                             d_i$sim.ElapsedTime <- d_i$sim.EndTime - d_i$sim.StartTime

                             saveRDS(d_i,
                                     file = here::here(sim.Path,
                                                       d_i$sim.File))

                           })

    })
    # Stop the cluster:
    snow::stopCluster(cl)

    return(t.snow)

  }


## @knitr do_fit_doFuture


do_fit_doFuture <-
  function(fit_refs,
           nClust = 48,
           nPROC = 1,
           save.directory = "self-sim",
           alternative.fit.Path = NULL,
           # model_what = "resid.random",
           clusterLOG.filename = paste0("fit_clusterLOG_",
                                        Sys.Date(),
                                        ".txt"),
           sleeptime = 1 # seconds to wait before running clusters
  ){

    source(here::here("functions",
                      "functions_Mplus.R"))



    debug <- TRUE

    d <- fit_refs

    registerDoFuture()

    plan("multisession")


    plyr::a_ply(d,
                # "uSeed",
                1,
                # base::transform,
                function(d_i){

                  # if(i<=nClust) Sys.sleep(sleeptime*i)

                  d_i <- as.list(d_i)

                  fit.StartTime <- Sys.time()

                  df <- readRDS(file = here::here(d_i$sim.Path,
                                                  d_i$sim.File))$sim.Dataset %>%
                    filter(subject <=d_i$N,
                           t <= d_i$T)
                  file.name <- gsub(".rds", "", d_i$fit.File)

                  print(paste(file.name,
                              "started at",
                              fit.StartTime)
                  )

                  tryRes <-
                    try(
                      output.fit <- run_MplusAutomation(df = df,
                                                        PROCESSORS = nPROC,
                                                        BITERATIONS.min = d_i$iter,
                                                        THIN = d_i$thin,
                                                        model_what = d_i$type,
                                                        file.name = file.name)
                    )

                  d_i$fit.Dataset <- tryRes # output.fit
                  d_i$fit.StartTime <- fit.StartTime
                  d_i$fit.EndTime <- Sys.time()
                  d_i$fit.ElapsedTime <- d_i$fit.EndTime - d_i$fit.StartTime

                  fit.Path <- ifelse(is.null(alternative.fit.Path),
                                     d_i$fit.Path,
                                     alternative.fit.Path)

                  saveRDS(d_i,
                          file = here::here(fit.Path,
                                            d_i$fit.File))

                  # return(paste("Running iteration:",
                  #              i,
                  #              " / ",
                  #              nrow(d),
                  #              "Time:",
                  #              as.character(Sys.time())))

                  print(paste(file.name,
                              "finished at",
                              d_i$fit.EndTime)
                  )                  # return(d_i)
                },
                .parallel = TRUE)


  }


## @knitr do_harvest_doFuture


do_harvest_doFuture <- function(fit.files){

  registerDoFuture()

  plan("multisession")

  results <- foreach(f = 1:length(fit.files),
                     .combine = rbind,
                     .errorhandling = 'remove') %dopar% {
                       fit_extract(fit.files[f])
                     }
  return(results)
}

## @knitr complete_pipe

source(here::here("functions",
                  "functions_self-sim-pipeline.R"))

#' We need to
#'  1. remake `sim_refs` to include NAR
#'  2. simulate NAR for N=100, T=100 in `simulation-files/sim-files`
#'  3. make reference table `fit_refs`
#'
#' Then need to subset the sim files for smaller N & T, and update the fit
#' reference table (`fit_refs`) such that the subset N & T values are included
#' in fit the file names (and respective Mplus files)
#'
#' I do this by:
#'  1. First making a `sim_refs` with smaller N & T (25, 50, 100) while
#' keeping uSeed intact, and then
#'  2. Change the `make_fit_refs` function such that the fit.File includes
#' N and T values, and finally
#'  3. Change `do_fit_parallel` to subset sim datasets prior to fitting.
#'
#' Note that we are no longer interested in DAR, so they should be left out.
#'


## Making the dataframe of simulation reference files

sim_refs <- make_sim_refs(conditions = list(T = c(30, 100),
                                                     N = c(100),
                                                     Model = c("BinAR",
                                                               "Chi2AR",
                                                               "DAR",
                                                               "PoDAR",
                                                               "NAR"),
                                                     l2.dist = c("Gaussian",
                                                                 "Chi2"),
                                                     phi = c(0.4)
),
save.directory = "simulation-files/sim-files",
Reps = 1000)

## Making a backup of the dataframe of simulation reference files

saveRDS(sim_refs,
        here::here("simulation-files/refs",
                   "sim-refs_Model-BinAR.ChiAR.DAR.PoDAR.NAR_N-100_T-30.100.rds")
)

## Reading the backup of the dataframe of simulation reference files

sim_refs.big <- readRDS(
  here::here(
    "simulation-files/refs",
    "sim-refs_Model-BinAR.ChiAR.DAR.PoDAR.NAR_N-100_T-30.100.rds"
  )
)

## Running the simulation

Sys.time()
system.time(
  t.sim <- do_sim_parallel(sim_refs = sim_refs_only.NAR,
                           save.directory = "simulation-files/sim-files"
  )
)
Sys.time()


## Making the dataframe of analysis output reference files

fit_refs <- make_fit_refs(sim_refs = sim_refs,
                          save.directory = "simulation-files/fit-files")


## Saving a backup of the dataframe of analysis output reference files

saveRDS(fit_refs,
        here::here("simulation-files/refs",
                   "fit-refs_Model-BinAR.ChiAR.PoDAR.NAR_N-25.50.100_T-25.50.100_iter-2000_thin-5_type-resid.fixed.random.rds")
)


## Reading the backup of the dataframe of analysis output reference files

fit_refs <- readRDS(here::here("simulation-files/refs",
                               "fit-refs_Model-BinAR.ChiAR.PoDAR.NAR_N-25.50.100_T-25.50.100_iter-2000_thin-5_type-resid.fixed.random.rds")
)


## Running the analysese in parallel

Sys.time()
system.time(
  t.fit <- do_fit_doFuture(fit_refs = fit_refs_remaining,
                           nClust = 48,
                           nPROC = 1,
                           sleeptime = 3,
                           save.directory = "simulation-files/fit-files"
  )
)
Sys.time()

## Harvesting the results in parallel

harvest.dir <- "simulation-files/fit-files"
# l.files <- list.files(path = here(harvest.dir),
#                       pattern = glob2rx("*_N-25_T-25*.rds"))
l.files <- list.files(path = here(harvest.dir),
                      pattern = glob2rx("*.rds"))


fit.files <- l.files %>%
  here(harvest.dir, .)


# registerDoFuture()
#
# plan("multisession")

Sys.time()
system.time(
  harv <- do_harvest_doFuture(fit.files)
)
Sys.time()


saveRDS(harv,
        "fit-harvest_64k.rds")
