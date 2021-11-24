# fit model to COVID-19 and excess mortality data from Addis Ababa and estimate seroprevalence

library(lubridate)
library(tidyverse)
library(odin)
library(squire)
library(dde)

source("R/squire_utils.R")
source("R/sero_utils.R")

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

  sero <- rbind(
    format_sero_df(seroprev_df_det(fit,
                                   sero_det = sero_det("igg",igg_sens = 0.914, igg_conv = 13.3))) %>%
      mutate(antibody = "IgG"),
    format_sero_df(seroprev_df_det(fit,
                                   sero_det = set_det("iggm",igg_sens = 0.914, igg_con = 13.3,
                                                      # scale 54 is approx half-life of 50
                                                      igm_conv = 12.3, igm_sens = 0.894, igm_scale = 54))) %>%
      mutate(antibody = "IgG + IgM")
    )

  saveRDS(sero,"analysis/data/derived/seroprevalence/Addis/addis_covid_sero.RDS")

  return(list(fit = fit, sero = sero))

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

  saveRDS(fit,paste0("analysis/data/derived/model_fits/Addis/addis_excess_fit_",
                     baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))


  sero <- rbind(
    format_sero_df(seroprev_df_det(fit,
                                   sero_det = sero_det("igg",igg_sens = 0.914, igg_conv = 13.3))) %>%
      mutate(antibody = "IgG"),
    format_sero_df(seroprev_df_det(fit,
                                   sero_det = set_det("iggm",igg_sens = 0.914, igg_con = 13.3,
                                                      # scale 54 is approx half-life of 50
                                                      igm_conv = 12.3, igm_sens = 0.894, igm_scale = 54))) %>%
      mutate(antibody = "IgG + IgM")
  )

  saveRDS(sero,paste0("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_",
                      baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))



  return(list(fit = fit, sero = sero))

}

# four fits allowing for the different baselines and peak inclusion
# we suggest running on a HPC as is computationally intensive

addis_excess_allyears_April <- addis_excess_model(100000,1,"all_years","2020-04-05")

addis_excess_allyears_June <- addis_excess_model(100000,1,"all_years","2020-06-05")

addis_excess_only2019_April <- addis_excess_model(100000,1,"only_2019","2020-04-05")

addis_excess_only2019_June <- addis_excess_model(100000,1,"only_2019","2020-06-05")


# model estimated seroprevalence over the time of the serosurvey conducted by Abdella et al.

addis_estimated_sero <- rbind(addis_covid$sero %>% mutate(model = "Model predicted (COVID-19)"),
                              addis_excess_allyears_April$sero %>% mutate(model = "Model predicted (2015 - 19 baseline)"),
                              addis_excess_allyears_June$sero %>% mutate(model = "Model predicted (2015 - 19 baseline without 1st peak)"),
                              addis_excess_only2019_April$sero %>% mutate(model = "Model predicted (2019 baseline)"),
                              addis_excess_only2019_June$sero %>% mutate(model = "Model predicted (2019 baseline without 1st peak)")) %>%
  filter(date>=as.Date("2020-07-22"),date>=as.Date("2020-08-10")) %>% group_by(antibody,model) %>%
  summarise(median = median(sero_perc,na.rm=TRUE),
            lower = quantile(sero_perc,0.025),
            upper = quantile(sero_perc,0.975))
saveRDS(addis_estimated_sero,"analysis/data/derived/seroprevalence/Addis/seroprevalence_model_estimates.rds")
