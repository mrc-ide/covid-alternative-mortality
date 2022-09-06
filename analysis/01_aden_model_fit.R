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


  sero_scenarios <- readRDS("analysis/data/derived/seroprevalence/Aden/sero_scenarios_without_IFR.RDS")

  sero_summary <- c()
  sero_ts <- c()

  for(i in 1:nrow(sero_scenarios)){

    sero_it <- seroprev_df_det(fit,
                               sero_det = sero_det(sero_scenarios$antibody_type[i],
                                                   igg_sens = 0.967,igm_sens = 1,
                                                   igg_scale = sero_scenarios$igg_serorev[i],
                                                   igm_scale = sero_scenarios$igm_serorev[i]))

    sero_summary_it <- sero_it %>% filter(date>=as.Date("2020-11-28"),
                                          date<=as.Date("2020-12-13")) %>%
      summarise(median = median(sero_perc),
                var = var(sero_perc)) %>%
      mutate(igg_scenario=sero_scenarios$igg_scenario[i],
             igm_scenario=sero_scenarios$igm_scenario[i],
             antibody=sero_scenarios$antibody_type[i],
             ifr=sero_scenarios$ifr[i])

    sero_summary <- rbind(sero_summary,sero_summary_it)


    sero_ts_out <- format_sero_df(sero_it) %>%
      mutate(igg_scenario=sero_scenarios$igg_scenario[i],
             igm_scenario=sero_scenarios$igm_scenario[i],
             antibody=sero_scenarios$antibody_type[i])

    sero_ts <- rbind(sero_ts,sero_ts_out)


  }

  return(list(fit = fit,sero_summary = sero_summary,sero_ts=sero_ts))
}

aden_covid <- aden_fit_covid(100000,1)

# suppress output so can upload to git
aden_covid$fit$output <- NULL
saveRDS(aden_covid$fit,"analysis/data/derived/model_fits/Aden/aden_covid_fit.rds")
saveRDS(aden_covid$sero_summary,"analysis/data/derived/seroprevalence/Aden/aden_sero_covid_summary.RDS")
saveRDS(aden_covid$sero_ts,"analysis/data/derived/seroprevalence/Aden/aden_sero_covid_time_series.RDS")


# fit to excess mortality - default fit
# no sero for this model as just to check IFR - sero is considered below when varying IFR

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


  return(list(fit = fit, fit_complete = fit_complete))

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

  ## now do sero

  sero_scenarios <- readRDS("analysis/data/derived/seroprevalence/Aden/sero_scenarios_without_IFR.RDS")

  sero_summary <- c()
  sero_ts <- c()
  sero_ll <- c()

  for(i in 1:nrow(sero_scenarios)){

       sero_it <- seroprev_df_det(fit_complete,
                                  sero_det = sero_det(sero_scenarios$antibody_type[i],
                                                      igg_sens = 0.967, igg_conv = 13.3,
                                                      igg_scale = sero_scenarios$igg_serorev[i],
                                                      igm_scale = sero_scenarios$igm_serorev[i],
                                                      igm_conv = 12.3, igm_sens = 1))

       sero_summary_it <- sero_it %>% filter(date>=as.Date("2020-11-28"),
                                             date<=as.Date("2020-12-13")) %>%
         summarise(median = median(sero_perc),
                   var = var(sero_perc)) %>%
         mutate(igg_scenario=sero_scenarios$igg_scenario[i],
                igm_scenario=sero_scenarios$igm_scenario[i],
                antibody=sero_scenarios$antibody_type[i],
                ifr=sero_scenarios$ifr[i])

       sero_summary <- rbind(sero_summary,sero_summary_it)


    sero_ts_out <- format_sero_df(sero_it) %>%
      mutate(igg_scenario=sero_scenarios$igg_scenario[i],
             igm_scenario=sero_scenarios$igm_scenario[i],
             antibody=sero_scenarios$antibody_type[i],
             ifr=IFR)

    sero_ts <- rbind(sero_ts,sero_ts_out)

    if(sero_scenarios$antibody_type[i]=="igg"){
      no_pos <- 500
    }else if(sero_scenarios$antibody_type[i]=="igm"){
      no_pos <- 4
    }else if(sero_scenarios$antibody_type[i]=="iggm"){
      no_pos <- 549
    }else{
      stop("No other antibody_type.")
    }

    sero_ll_out <- data.frame("ll"=sero_log_likelihood_aden(sero_it,
                                                            obs_x=no_pos,obs_n=2001)) %>%
      mutate(igg_scenario=sero_scenarios$igg_scenario[i],
             igm_scenario=sero_scenarios$igm_scenario[i],
             antibody=sero_scenarios$antibody_type[i],
             ifr=IFR)

    sero_ll <- rbind(sero_ll,sero_ll_out)


  }

  return(list(fit_complete = fit_complete, sero_summary = sero_summary, sero_ts = sero_ts, sero_ll = sero_ll))

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


## combine sero summaries
excess_sero_summary <- rbind(aden_excess_IFR02$sero_summary,
                             aden_excess_IFR03$sero_summary,
                             aden_excess_IFR04$sero_summary,
                             aden_excess_IFR05$sero_summary)
saveRDS(excess_sero_ts,"analysis/data/dervied/seroprevalence/Aden/aden_sero_excess_summary.RDS")



## combine sero time series
excess_sero_ts <- rbind(aden_excess_IFR02$sero_ts,
                        aden_excess_IFR03$sero_ts,
                        aden_excess_IFR04$sero_ts,
                        aden_excess_IFR05$sero_ts)
saveRDS(excess_sero_ts,"analysis/data/dervied/seroprevalence/Aden/aden_sero_excess_time_series.RDS")


## combine sero log likelihood
excess_sero_ll <- rbind(aden_excess_IFR02$sero_ll,
                        aden_excess_IFR03$sero_ll,
                        aden_excess_IFR04$sero_ll,
                        aden_excess_IFR05$sero_ll)
saveRDS(excess_sero_ll,"analysis/data/dervied/seroprevalence/Aden/aden_sero_excess_ll.RDS")





