# fit model to mortality data from Aden and estimate seroprevalence

library(tidyverse)
library(squire)

source("R/squire_utils.R")
source("R/sero_utils.R")

# fit to COVID-19 deaths

aden_fit_covid <- function(n_mcmc,rf){

  aden_deaths <- readRDS("analysis/data/derived/deaths_time_series/aden_covid_deaths.RDS")

  aden_pop <- readRDS("analysis/data/raw/aden_population.RDS")

  fit <- fit_spline_rt(data=aden_deaths,
                       country="Yemen",
                       population=round(aden_pop$prop*1000000),
                       reporting_fraction=rf,
                       n_mcmc = n_mcmc,
                       replicates = 100,
                       rw_duration = 14,
                       hosp_beds = 1133,
                       icu_beds = 17)

  #saveRDS(fit,"analysis/data/derived/model_fits/Aden/aden_covid_fit.rds")

  fit_sero <- rbind(
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 109.5,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="100 days"),
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 310,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="280 days")
  )
  saveRDS(fit_complete_sero,"analysis/data/derived/seroprevalence/Aden/aden_covid_sero.rds")

  return(list(fit = fit, fit_sero = fit_sero))
}

aden_covid <- aden_fit_covid(100000,1)

# suppress output so can upload to git
aden_covid$fit$output <- NULL
saveRDS(aden_covid$fit,"analysis/data/derived/model_fits/Aden/aden_covid_fit.rds")


# fit to excess mortality


aden_fit_excess <- function(n_mcmc,rf){

  aden_deaths <- readRDS("analysis/data/derived/deaths_time_series/aden_excess_deaths.rds") %>%
    mutate(deaths = as.integer(ifelse(deaths<0,0,deaths)))

  aden_pop <- readRDS("analysis/data/raw/aden_population.RDS")

  fit <- fit_spline_rt(data=aden_deaths,
                       country="Yemen",
                       population=round(aden_pop$prop*1000000),
                       reporting_fraction=rf,
                       n_mcmc = n_mcmc,
                       replicates = 100,
                       rw_duration = 14,
                       hosp_beds = 1133,
                       icu_beds = 17)

  #saveRDS(fit,"analysis/data/derived/model_fits/Aden/aden_excess_fit_raw.rds")

  # then project forward to the end of 2020 assuming Rt stays as estimated at end of fit
  fit_complete <- squire::projections(fit,time_period = 125)
  #saveRDS(fit_complete,"analysis/data/derived/model_fits/Aden/aden_excess_fit_complete.rds")


  # seroprevalence under different assumptions of IgG sero-reversion half-lives

  fit_complete_sero <- rbind(
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 155,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="140 days"),
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 176.5,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="160 days"),
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 199,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="180 days"),
    format_sero_df(seroprev_df_det(fit_complete,
                                   sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,igg_scale = 219.67,
                                                       igm_conv = 12.3, igm_sens = 1)))
    %>% mutate(scenario="200 days")
  )
  saveRDS(fit_complete_sero,"analysis/data/derived/seroprevalence/Aden/aden_excess_complete_sero.rds")


  return(list(fit = fit, fit_complete = fit_complete, fit_complete_sero = fit_complete_sero))

}

# fit with 100000 mcmc iterations
# we suggest running on a HPC as is computationally intensive

aden_excess <- aden_fit_excess(100000,1)

# suppress output so can upload to git
aden_excess$fit$output <- NULL
saveRDS(fit,"analysis/data/derived/model_fits/Aden/aden_excess_fit_raw.rds")

# vary IFR in model fit

aden_fit_excess_vary_IFR <- function(n_mcmc,rf,prob_death_scale,IFR){

  aden_deaths <- readRDS("analysis/data/derived/deaths_time_series/aden_excess_deaths.rds") %>%
    mutate(deaths = as.integer(ifelse(deaths<0,0,deaths)))

  aden_pop <- readRDS("analysis/data/raw/aden_population.RDS")

  fit <- fit_spline_rt_ifr_change(data=aden_deaths,
                                  country="Yemen",
                                  population=round(aden_pop$prop*1000000),
                                  reporting_fraction=rf,
                                  n_mcmc = n_mcmc,
                                  replicates = 100,
                                  rw_duration = 14,
                                  hosp_beds = 1e10,
                                  icu_beds = 1e10,
                                  prob_death_scale = prob_death_scale)

  fit_complete <- squire::projections(fit,time_period = 125)
  #saveRDS(fit_complete,paste0("analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR0",IFR*10,".rds"))

  # log likelihood of seroprevalence under different assumptions of IgG sero-reversion half-life

  ll <- rbind(
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                        sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                            igg_scale = 109.5,
                                                            igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=100,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 132,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=120,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 155,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=140,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 176.5,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=160,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 199,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=180,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 219.67,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=200,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 239,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=220,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 266.67,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=240,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 287.5,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=260,"ifr"=IFR),
    c("ll"=sero_log_likelihood_aden(seroprev_df_det(fit_complete,
                                               sero_det = sero_det("iggm", igg_sens = 0.967, igg_conv = 13.3,
                                                                   igg_scale = 310,
                                                                   igm_conv = 12.3, igm_sens = 1)),obs_x=549,obs_n=2001),
      "half_life"=280,"ifr"=IFR)
  )



  saveRDS(ll,paste0("analysis/data/derived/seroprevalence/Aden/aden_ll_IFR0",IFR*10,".rds"))

  return(list(fit_complete = fit_complete, ll = ll))



}

aden_excess_IFR02 <- aden_fit_excess_vary_IFR(100000,1,0.2/0.2,0.2)
aden_excess_IFR03 <- aden_fit_excess_vary_IFR(100000,1,0.3/0.2,0.3)
aden_excess_IFR04 <- aden_fit_excess_vary_IFR(100000,1,0.4/0.2,0.4)
aden_excess_IFR05 <- aden_fit_excess_vary_IFR(100000,1,0.6/0.2,0.5)

# suppress output so can upload to git
aden_excess_IFR02$fit_complete$output <- NULL
aden_excess_IFR03$fit_complete$output <- NULL
aden_excess_IFR04$fit_complete$output <- NULL
aden_excess_IFR05$fit_complete$output <- NULL

saveRDS(aden_excess_IFR02$fit_complete,
        "analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR02.rds")
saveRDS(aden_excess_IFR03$fit_complete,
        "analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR03.rds")
saveRDS(aden_excess_IFR04$fit_complete,
        "analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR04.rds")
saveRDS(aden_excess_IFR05$fit_complete,
        "analysis/data/derived/model_fits/Aden/aden_excess_fit_complete_IFR05.rds")








