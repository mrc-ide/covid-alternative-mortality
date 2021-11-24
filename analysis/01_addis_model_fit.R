# fit model to COVID-19 and excess mortality data from Addis Ababa

library(lubridate)
library(tidyverse)
library(odin)
library(squire)
library(dde)

source("R/squire_utils.R")
source("R/")

### fit to COVID-19 deaths

addis_covid_model <- function(n_mcmc,rf){

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

  saveRDS(fit,"analysis/data/derived/model_fits/Addis/addis_covid_fit.RDS")





  return(fit)

}

# fit with 100000 mcmc iterations
# we suggest running on a HPC as is computationally intensive

addis_covid <- addis_covid_model(100000,1)


### fit to excess mortality with different baselines and with and without May 2020 peak

addis_excess_model <- function(n_mcmc,rf,baseline_select,date_filter){

  addis_deaths <- readRDS("analysis/data/derived/addis_excess_deaths.RDS") %>%
    filter(baseline==baseline_select) %>%  select(dateburialgc,weekly) %>% rename(date=dateburialgc,deaths=weekly) %>%
    mutate(deaths=round(deaths),date=as.Date(date,format="%Y-%m-%d")) %>%
    filter(date >= as.Date(date_filter)&date <= as.Date("2021-01-01")) %>%
    mutate(deaths=ifelse(deaths<0,0,deaths))

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

  saveRDS(fit,paste0("analysis/data/derived/model_fits/Addisaddis_excess_fit_",
                     baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))

  return(fit)

}

# four fits allowing for the different baselines and peak inclusion
# we suggest running on a HPC as is computationally intensive

addis_excess_fit_allyears_April <- addis_excess_model(100000,1,"all_years","2020-04-05")

addis_excess_fit_allyears_June <- addis_excess_model(100000,1,"all_years","2020-06-05")

addis_excess_fit_only2019_April <- addis_excess_model(100000,1,"only_2019","2020-04-05")

addis_excess_fit_only2019_June <- addis_excess_model(100000,1,"only_2019","2020-06-05")



