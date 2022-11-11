## @knitr load_packages

libs_used <-
  c("tidyverse",
    "knitr",
    "MplusAutomation",
    "doSNOW",
    "doFuture",
    "here",
    "latex2exp",
    "cowplot",
    "glue",
    "patchwork",
    "ggExtra",
    "primes")

libs_needed <- libs_used[!libs_used %in% installed.packages()]

sapply(libs_needed, install.packages, dependencies = TRUE)
sapply(libs_used, require, character = TRUE)
