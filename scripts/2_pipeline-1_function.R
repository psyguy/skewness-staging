## @knitr make_sim_refs

make_sim_refs <-
  function(conditions = list(T = c(30, 100),
                             N = c(100),
                             Model = c("BinAR", "Chi2AR", "DAR"),
                             l2.dist = c("Gaussian", "Chi2"),
                             phi = c(0.4)
  ),
  Reps = 100,
  simSeed = 0,
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

    conditions$simSeed <- simSeed
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

    # getting rid of factors and making them into strings
    factor.columns <- sapply(d.headers, is.factor)
    d.headers[factor.columns] <-
      sapply(d[factor.columns], as.character)

    # making uSeed numeric
    d.headers$uSeed <- d.headers$uSeed %>% as.numeric()

    # Save the references dataframe to a file, if desired
    if(!is.null(save.refs.filename)){
      saveRDS(d.headers,
              file = here::here(save.directory, paste0(save.refs.filename,
                                                       ".rds"))
      )
      write.csv(d.headers,
                file = here::here(save.directory, paste0(save.refs.filename,
                                                         ".csv")),
                row.names = FALSE)
    }

    # return the reference dataframe
    return(d.headers)
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
           sleeptime = 1) {
    # cl <- snow::makeSOCKcluster(nClust,
    #                             outfile = here::here(save.directory,
    #                                                  clusterLOG.filename))
    debug <- TRUE

    d <- sim_refs

    registerDoFuture()

    plan("multisession")

    ## To make sure the required functions are loaded on each cluster
    library(tidyverse)
    source(here::here("scripts",
                      "1_component-1_simulation.R"))
    plyr::a_ply(d,
                1,
                function(d_i) {
                  # if (i <= nClust)
                  #   Sys.sleep(sleeptime * i)

                  d_i <- as.list(d_i)

                  arguments <-
                    as.list(d_i[2:(length(d_i) - 4)])
                  arguments$seed <- d_i$uSeed

                  sim.Path <-
                    ifelse(is.null(alternative.sim.Path),
                           d_i$sim.Path,
                           alternative.sim.Path)

                  sim.StartTime <- Sys.time()

                  output_make_datasets <-
                    do.call(make_datasets,
                            arguments)

                  if (is.data.frame(output_make_datasets)) {
                    d_i$sim.Dataset <- output_make_datasets
                  } else{
                    # Adapting for v2 implementation
                    d_i$sim.Dataset <-
                      output_make_datasets$dataset
                    d_i$given.Means <-
                      output_make_datasets$Means
                    d_i$given.Phis <-
                      output_make_datasets$Phis
                  }
                  d_i$sim.StartTime <- sim.StartTime
                  d_i$sim.EndTime <- Sys.time()
                  d_i$sim.ElapsedTime <-
                    d_i$sim.EndTime - d_i$sim.StartTime

                  saveRDS(d_i,
                          file = here::here(sim.Path,
                                            d_i$sim.File))

                })

  }


## @knitr do_fit_doFuture

do_fit_doFuture <- function(fit_refs,
                            nClust = 48,
                            nPROC = 1,
                            save.directory = "self-sim",
                            alternative.fit.Path = NULL,
                            # model_what = "resid.random",
                            clusterLOG.filename = paste0("fit_clusterLOG_",
                                                         Sys.Date(),
                                                         ".txt"),
                            sleeptime = 1,
                            version = "v2") {

  ## To make sure the required functions are loaded on each cluster
  library(tidyverse)
  source(here::here("scripts",
                    "1_component-2_analysis.R"))




  debug <- TRUE

  d <- fit_refs

  registerDoFuture()

  plan("multisession")


  plyr::a_ply(d,
              1,
              function(d_i) {
                # if (i <= nClust)
                #   Sys.sleep(sleeptime * i)

                d_i <- as.list(d_i)

                fit.StartTime <-
                  Sys.time()

                sim.file <- readRDS(file = here::here(d_i$sim.Path,
                                                      d_i$sim.File))

                df <- sim.file$sim.Dataset %>%
                  filter(subject <= d_i$N,
                         t <= d_i$T)

                # Adapting for v2 implementation
                df_sum <- df %>%
                  group_by(subject) %>%
                  summarise(sample.Means = mean(x),
                            sample.Variances = var(x),
                            sample.Skewnesses = moments::skewness(x))

                file.name <-
                  gsub(".rds", "", d_i$fit.File)

                print(paste(file.name,
                            "started at",
                            fit.StartTime))

                tryRes <-
                  try(output.fit <-
                        run_MplusAutomation(
                          df = df,
                          PROCESSORS = nPROC,
                          BITERATIONS.min = d_i$iter,
                          THIN = d_i$thin,
                          model_what = d_i$type,
                          file.name = file.name
                        ))

                d_i$fit.Dataset <-
                  tryRes

                # Adapting for v2 implementation
                d_i$given.Means <- sim.file$given.Means %>% head(d_i$N)
                d_i$given.Phis <- sim.file$given.Phis %>% head(d_i$N)

                # Adapting for v2 implementation
                d_i$sample.Means <- df_sum$sample.Means
                d_i$sample.Variances <- df_sum$sample.Variances
                d_i$sample.Skewnesses <- df_sum$sample.Skewnesses

                d_i$fit.StartTime <-
                  fit.StartTime
                d_i$fit.EndTime <-
                  Sys.time()
                d_i$fit.ElapsedTime <-
                  d_i$fit.EndTime - d_i$fit.StartTime

                fit.Path <-
                  ifelse(is.null(alternative.fit.Path),
                         d_i$fit.Path,
                         alternative.fit.Path)

                saveRDS(d_i,
                        file = here::here(fit.Path,
                                          d_i$fit.File))

                print(paste(file.name,
                            "finished at",
                            d_i$fit.EndTime))
              },
              .parallel = TRUE)


}


## @knitr do_harvest_doFuture


do_harvest_doFuture <- function(fit.files,
                                nClust = 46,
                                sleeptime = 1){

  registerDoFuture()

  plan("multisession")

  results <- foreach(i = 1:length(fit.files),
                     .combine = rbind,
                     .errorhandling = 'remove') %dopar% {
                       if (i <= nClust)
                         Sys.sleep(sleeptime * i)
                       fit_extract(fit.files[i])
                     }
  return(results)
}


## @knitr pipeline_make_references

# sim_refs_base <- make_sim_refs(
#   conditions = list(
#     T = c(30, 100),
#     N = c(100),
#     Model = c("BinAR",
#               "Chi2AR",
#               "DAR",
#               "PoDAR",
#               "NAR"),
#     l2.dist = c("Gaussian",
#                 "Chi2"),
#     phi = c(0.4)
#   ),
#   simSeed = 0,
#   Reps = 1000,
#   save.directory = dir_simulation
# ) %>%
#   filter(T == 100, Model != "DAR")
#
# ## Adding another values for N and T; see the text
#
# sim_refs <- NULL
#
# for (TT in c(25, 50, 100)) {
#   for (NN in c(25, 50, 100)) {
#     sim_refs <- sim_refs_base %>%
#       mutate(N = NN,
#              T = TT) %>%
#       rbind(sim_refs)
#   }
# }
#
# ## Saving a backup of the dataframe of references of simulation files
#
# saveRDS(sim_refs,
#         here::here(dir_references,
#                    "sim-refs.rds"))
#
# ## Making the dataframe of analysis output reference files
#
# fit_refs <- make_fit_refs(sim_refs = sim_refs,
#                           save.directory = dir_analysis)
#
#
# ## Saving a backup of the dataframe of references of analysis output files
#
# saveRDS(fit_refs,
#         here::here(dir_references,
#                    "fit-refs.rds"))


## @knitr pipeline_set_directories

## Root directory of simulation study files
dir_files <- "simulation-files"

## Where reference tables are saved
dir_references <- paste(dir_files,
                        "refs",
                        sep = "/")

## Where simulated datasets are saved
dir_simulation <- paste(dir_files,
                        "sim-files",
                        sep = "/")

## Where analysis output files are saved
dir_analysis <- paste(dir_files,
                      "fit-files",
                      sep = "/")

## Where harvested study results are saved
dir_harvest <- paste(dir_files,
                     "harvest-files",
                     sep = "/")

## where plots are saved

dir_plots <- "figures"
