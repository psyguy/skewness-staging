## @knitr run_MplusAutomation

run_MplusAutomation <- function(df,
                                PROCESSORS = 2,
                                CHAINS = 2,
                                THIN = 2,
                                BITERATIONS.min = 2000,
                                BITERATIONS.max =
                                  BITERATIONS.min * BITERATION.minmax.factor,
                                BITERATION.minmax.factor = 2.5,
                                out.folder = "Mplus-files/",
                                model_what = "resid.random",
                                file.name = Sys.Date()) {
  inp.name <- paste0(out.folder,
                     file.name)

  MplusAutomation::prepareMplusData(df,
                                    paste0(inp.name, ".dat"))

  VARIABLE <- glue::glue("CLUSTER = subject;
                      LAGGED = x(1);
                      TINTERVAL = t(1);")

  ANALYSIS  <- glue::glue(
    "TYPE = TWOLEVEL RANDOM;
  	ESTIMATOR = BAYES;
  	PROCESSORS = {PROCESSORS};
    CHAINS = {CHAINS};
    THIN = {THIN};
  	BITERATIONS = {BITERATIONS.max}({BITERATIONS.min});"
  )

  PLOT <- glue::glue("TYPE = PLOT3;
                     FACTORS = ALL (500);")


  if (model_what == "resid.random")
    model_string <- glue::glue("%WITHIN%
  	phi | x ON x&1;
  	logv | x;
  	%BETWEEN%
  	x phi logv WITH x phi logv;")

  if (model_what == "resid.fixed")
    model_string <- glue::glue("%WITHIN%
  	phi | x ON x&1;
  	%BETWEEN%
  	x phi WITH x phi;")


  if (model_what == "within.between") {
    model_string <- glue::glue("%WITHIN%
  	x;
  	%BETWEEN%
  	x;")
    VARIABLE <- "CLUSTER = subject;"
    PLOT <- "TYPE = PLOT3;"

  }


  model.ar1 <- MplusAutomation::mplusObject(
    TITLE = inp.name,
    rdata = df,
    usevariables = c("subject", "t", "x"),
    VARIABLE = VARIABLE,
    ANALYSIS = ANALYSIS,
    MODEL = model_string,
    OUTPUT = "TECH1 TECH2 TECH3 TECH8 FSCOMPARISON STANDARDIZED STDYX STDY;",
    PLOT = PLOT
  )


  fit.ar1 <- MplusAutomation::mplusModeler(
    model.ar1,
    check = FALSE,
    modelout = paste0(inp.name, ".inp"),
    hashfilename = FALSE,
    run = 1L
  )

  return(fit.ar1)
}
