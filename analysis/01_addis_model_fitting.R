# fit model to COVID-19 and excess mortality from Addis Ababa

library(lubridate)
library(tidyverse)
library(odin)
library(squire)
library(dde)

source("R/squire_utils.R")

### COVID-19 ###

addis_covid_fit <- function(n_mcmc,rf){

  addis_deaths <- readRDS("analysis/data/derived/addis_covid_deaths.RDS") %>%
    filter(date >= as.Date("05/04/2020"))

  addis_pop <- readRDS("analysis/data/raw/addis_population.RDS")

  fit <- fit_spline_rt(data=addis_deaths,
                       country="Ethiopia",
                       population=round(addis_pop$prop*4800000),
                       reporting_fraction=rf,
                       n_mcmc = n_mcmc,
                       replicates = 100,
                       rw_duration = 14,
                       hosp_beds = 6867,
                       icu_beds = 103)

  saveRDS(fit,"analysis/data/derived/addis_covid_fit.RDS")

  ### do we want seroprevalence in this function

  seroprev <- seroprev_df(fit)

  saveRDS(seroprev,paste0("N:/Ruth/covid-global-mortality/scripts/unfinished_scripts/Addis/fits_oct25/seroprev/sero_addis_covid_",n_mcmc,"_rf",rf,"_",Sys.Date(),".RDS"))

  return(list(fit=fit,seroprev=seroprev))

}






