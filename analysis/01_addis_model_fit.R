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

  addis_deaths <- readRDS("analysis/data/derived/deaths_time_series/addis_covid_deaths.RDS") %>%
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

  #saveRDS(fit,"analysis/data/derived/model_fits/Addis/addis_covid_fit.RDS")

  sero_igg <- seroprev_df_det(fit,
                              sero_det = sero_det("igg", igg_sens = 0.914, igg_conv = 13.3))

  sero_iggm <- seroprev_df_det(fit,
                               sero_det = sero_det("iggm", igg_sens = 0.914, igg_conv = 13.3,
                                                   igm_conv = 12.3, igm_sens = 0.894, igm_scale = 54))

  sero <- rbind(format_sero_df(sero_igg) %>% mutate(antibody = "IgG"),
                format_sero_df(sero_iggm) %>% mutate(antibody = "Combined IgG/IgM"))

  saveRDS(sero,"analysis/data/derived/seroprevalence/Addis/addis_covid_sero.RDS")


  sero_igg_summary <- sero_igg %>%
    filter(date>=as.Date("2020-07-22"),
           date<=as.Date("2020-08-10")) %>%
    summarise(median = median(sero_perc),
              var = var(sero_perc))

  sero_iggm_summary <- sero_iggm %>%
    filter(date>=as.Date("2020-07-22"),
           date<=as.Date("2020-08-10")) %>%
    summarise(median = median(sero_perc),
              var = var(sero_perc))


  sero_model_est <- rbind(sero_igg_summary %>% mutate(antibody=="IgG"),
                          sero_iggm_summary %>% mutate(antibody=="Combined IgG/IgM"))


  #saveRDS(sero_model_est,"analysis/data/derived/seroprevalence/Addis/addis_covid_sero_model_est.RDS")

  return(list(fit = fit, sero = sero, sero_model_est = sero_model_est))

}

# fit with 100000 mcmc iterations
# we suggest running on a HPC as is computationally intensive

addis_covid <- addis_covid_model(100000,1)

# suppress output so can upload to git
addis_covid$fit$output <- NULL
saveRDS(addis_covid$fit,"analysis/data/derived/model_fits/Addis/addis_covid_fit.RDS")


### fit to excess mortality with different baselines and with and without May 2020 peak

addis_excess_model <- function(n_mcmc,rf,baseline_select,date_filter){

  addis_deaths <- readRDS("analysis/data/derived/deaths_time_series/addis_excess_deaths.RDS") %>%
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

  # saveRDS(fit,paste0("analysis/data/derived/model_fits/Addis/addis_excess_fit_",
  #                    baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))


  sero_igg <- seroprev_df_det(fit,
                              sero_det = sero_det("igg", igg_sens = 0.914, igg_conv = 13.3))

  sero_iggm <- seroprev_df_det(fit,
                               sero_det = sero_det("iggm", igg_sens = 0.914, igg_conv = 13.3,
                                                   igm_conv = 12.3, igm_sens = 0.894, igm_scale = 54))

  sero <- rbind(format_sero_df(sero_igg) %>% mutate(antibody = "IgG"),
                format_sero_df(sero_iggm) %>% mutate(antibody = "Combined IgG/IgM"))

  saveRDS(sero,paste0("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_",
                      baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))


  sero_igg_summary <- sero_igg %>%
    filter(date>=as.Date("2020-07-22"),
           date<=as.Date("2020-08-10")) %>%
    summarise(median = median(sero_perc),
              var = var(sero_perc))

  sero_iggm_summary <- sero_iggm %>%
    filter(date>=as.Date("2020-07-22"),
           date<=as.Date("2020-08-10")) %>%
    summarise(median = median(sero_perc),
              var = var(sero_perc))


  sero_model_est <- rbind(sero_igg_summary %>% mutate(antibody=="IgG"),
                          sero_iggm_summary %>% mutate(antibody=="Combined IgG/IgM"))


  #
  # sero_model_est <- rbind(sero_igg %>% mutate(antibody = "IgG"),
  #                         sero_iggm %>% mutate(antibody = "Combined IgG/IgM")) %>%
  #   filter(date>=as.Date("2020-07-22")&date<=as.Date("2020-08-10")) %>% group_by(antibody) %>%
  #   summarise(median=median(sero_perc,na.rm=TRUE),
  #             lower=quantile(sero_perc,0.025),
  #             upper=quantile(sero_perc,0.975))
  #
  # saveRDS(sero_model_est,paste0("analysis/data/derived/seroprevalence/Addis/addis_excess_sero_model_est_",
  #                               baseline_select,"_",format.Date(date_filter,"%B"),".RDS"))

  return(list(fit = fit, sero = sero, sero_model_est = sero_model_est))

}

# four fits allowing for the different baselines and peak inclusion
# we suggest running on a HPC as is computationally intensive

addis_excess_all_years_April <- addis_excess_model(100000,1,"all_years","2020-04-05")

addis_excess_all_years_June <- addis_excess_model(100000,1,"all_years","2020-06-05")

addis_excess_only_2019_April <- addis_excess_model(100000,1,"only_2019","2020-04-05")

addis_excess_only_2019_June <- addis_excess_model(100000,1,"only_2019","2020-06-05")


# model estimated seroprevalence over the time of the serosurvey conducted by Abdella et al.

addis_sero_model_est <- rbind(addis_covid$sero_model_est %>% mutate(model = "Model predicted (COVID-19)"),
                              addis_excess_all_years_April$sero_model_est %>%
                                mutate(model = "Model predicted (2015 - 19 baseline)"),
                              addis_excess_all_years_June$sero_model_est %>%
                                mutate(model = "Model predicted (2015 - 19 baseline without 1st peak)"),
                              addis_excess_only_2019_April$sero_model_est %>%
                                mutate(model = "Model predicted (2019 baseline)"),
                              addis_excess_only_2019_June$sero_model_est %>%
                                mutate(model = "Model predicted (2019 baseline without 1st peak)"))

saveRDS(addis_sero_model_est,"analysis/data/derived/seroprevalence/Addis/addis_seroprevalence_model_estimates.rds")

# suppress output so can upload to git
addis_excess_all_years_April$fit$output <- NULL
addis_excess_all_years_June$fit$output <- NULL
addis_excess_only_2019_April$fit$output <- NULL
addis_excess_only_2019_June$fit$output <- NULL

saveRDS(addis_excess_all_years_April$fit,
        "analysis/data/derived/model_fits/Addis/addis_excess_fit_all_years_April.RDS")
saveRDS(addis_excess_all_years_June$fit,
        "analysis/data/derived/model_fits/Addis/addis_excess_fit_all_years_June.RDS")
saveRDS(addis_excess_only_2019_April$fit,
        "analysis/data/derived/model_fits/Addis/addis_excess_fit_only_2019_April.RDS")
saveRDS(addis_excess_only_2019_June$fit,
        "analysis/data/derived/model_fits/Addis/addis_excess_fit_only_2019_June.RDS")



