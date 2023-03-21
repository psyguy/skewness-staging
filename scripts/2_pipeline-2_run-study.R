
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
  filter(!(sim.File %in% sim_refs_done)) %>%
  filter(T == 100, Model != "DAR")
# %>%
  # filter(Rep <= 500)

## Simulating the remaining datasets

Sys.time()
system.time(
  t.sim <- do_sim_parallel(sim_refs = sim_refs_remaining,
                           nClust = 46,
                           save.directory = dir_simulation)
)
Sys.time()


# Analysis part -----------------------------------------------------------

## Reading the backup of the dataframe of references of analysis output files

fit_refs <- readRDS(here::here(dir_references,
                               "fit-refs.rds"))
# %>%
  # filter(T == 100, Model != "DAR") %>%
  # filter(Rep > 500)

## Making a list of analysis output files

length(fit_files_done <- list.files(here::here(dir_analysis),
                             pattern = "fit_uSeed-"))

## Finding the analyses that are yet to be done

fit_refs_remaining <- fit_refs %>%
  filter(!(fit.File %in% fit_files_done))

## Running the analyses of remaining datasets in parallel

Sys.time()
system.time(
  t.fit <- do_fit_doFuture(
    fit_refs = fit_refs_remaining,
    nClust = 46,
    nPROC = 1,
    sleeptime = 3,
    save.directory = dir_analysis
  )
)
Sys.time()

fit_refs %>%
  filter(Rep %in% c(1:10)) %>%
  pull(fit.File) %>%
  here(dir_analysis, .) %>%
  zip::zipr("Rep_1-10.zip", .)

chunk_size <- 10
for (chunk in 1:(1000/chunk_size)) {

  index_range <- (chunk-1)*chunk_size + (1:chunk_size)
  print(ts <- Sys.time())

  print(zip_name <- paste0("Rep_",
                           index_range %>% head(1),
                           "-",
                           index_range %>% tail(1),
                           ".zip"))

  sims <- fit_refs %>%
    filter(Rep %in% index_range) %>%
    pull(sim.File) %>%
    here(dir_simulation, .)

  zip::zipr(paste0("archive_sim_",
                   zip_name) %>%
              here(dir_files,
                   "archived-files",
                   .),
            sims)

  fits <- fit_refs %>%
    filter(Rep %in% index_range) %>%
    pull(fit.File) %>%
    here(dir_analysis, .)

  zip::zipr(paste0("archive_fit_",
                   zip_name) %>%
              here(dir_files,
                   "archived-files",
                   .),
            fits)

  te <- Sys.time()

  print(te - ts)

  }

## @knitr pipeline_harvest_results

source("scripts/0_requirements.R")
source("scripts/1_component-1_simulation.R")
source("scripts/1_component-2_analysis.R")
source("scripts/1_component-3_harvesting.R")
source("scripts/1_component-4_reporting.R")
source("scripts/2_pipeline-1_function.R")

## Harvesting the results in parallel

l.files <- list.files(path = here(dir_analysis),
                      pattern = glob2rx("*.rds"))


fit.files <- here::here(dir_analysis, l.files)


# list.files(here::here(dir_harvest),
#            pattern = "harvest-raw_") %>%
#   sort(decreasing = TRUE) %>%
#   paste(here::here(dir_harvest),
#         .,
#         sep = "/") %>%
#   readRDS()

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
