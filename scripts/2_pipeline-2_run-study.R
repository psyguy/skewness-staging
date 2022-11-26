
## @knitr pipeline_run_study

# Simulation part ---------------------------------------------------------

## Reading the backup of the dataframe of references of simulation files

sim_refs <- readRDS(here::here(dir_references,
                               "sim-refs.rds"))

## Making a list of already simulated datasets

sim_refs_done <- list.files(here::here(dir_simulation),
                            pattern = "sim_uSeed-")

## Finding the datasets that are yet to be simulated

sim_refs_remaining <- sim_refs %>%
  filter(!(sim.File %in% sim_refs_done))

## Simulating the remaining datasets

Sys.time()
system.time(
  t.sim <- do_sim_parallel(sim_refs = sim_refs_remaining,
                           save.directory = dir_simulation)
)
Sys.time()


# Analysis part -----------------------------------------------------------

## Reading the backup of the dataframe of references of analysis output files

fit_refs <- readRDS(here::here(dir_references,
                               "fit-refs.rds"))

## Making a list of analysis output files

fit_files_done <- list.files(here::here(dir_simulation),
                             pattern = "fit_uSeed-")

## Finding the analyses that are yet to be done

fit_refs_remaining <- fit_refs %>%
  filter(!(fit.File %in% fit_files_done))

## Running the analyses of remaining datasets in parallel

Sys.time()
system.time(
  t.fit <- do_fit_doFuture(
    fit_refs = fit_refs_remaining,
    nClust = 48,
    nPROC = 1,
    sleeptime = 3,
    save.directory = dir_analysis
  )
)
Sys.time()


## @knitr pipeline_harvest_results

## Harvesting the results in parallel

l.files <- list.files(path = here(dir_analysis),
                      pattern = glob2rx("*.rds"))


fit.files <- l.files %>%
  paste(here::here(dir_analysis),
        .,
        sep = "/")


list.files(here::here(dir_harvest),
           pattern = "harvest-raw_") %>%
  sort(decreasing = TRUE) %>%
  paste(here::here(dir_harvest),
        .,
        sep = "/") %>%
  readRDS()

Sys.time()
system.time(
  harvest_raw <- do_harvest_doFuture(fit.files)
)
Sys.time()


saveRDS(harvest_raw,
        here::here(dir_harvest,
                   paste0(
                     "harvest-raw_",
                     format(Sys.time(),
                            "%Y-%m-%d_%H-%M"),
                     ".rds"
                   )))


d_abridged <- harvest_raw %>%
  harvest_cleanup(return.abridged = TRUE)

saveRDS(d_abridged ,
        here::here(dir_harvest,
                   paste0(
                     "harvest-abridged_",
                     format(Sys.time(),
                            "%Y-%m-%d_%H-%M"),
                     ".rds"
                   )))


d_important <- harvest_raw %>% harvest_cleanup()

saveRDS(d_important ,
        here::here(dir_harvest,
                   paste0(
                     "harvest-important_",
                     format(Sys.time(),
                            "%Y-%m-%d_%H-%M"),
                     ".rds"
                   )))
