# fit model to excess mortality data from Aden and estimate seroprevalence

library(tidyverse)
library(squire)

source("R/squire_utils.R")
source("R/sero_utils.R")

aden_fit_excess_smooth <- function(n_mcmc,rf){

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

  saveRDS(fit,"analysis/data/derived/model_fits/Aden/aden_excess_fit_raw.RDS")

  # then project forward to the end of 2020 assuming Rt stays as estimated at end of fit
  fit_complete <- squire::projections(fit,time_period = 125)
  saveRDS(fit_complete,"analysis/data/derived/model_fits/Aden/aden_excess_fit_complete.RDS")


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
  saveRDS(fit_complete_sero,"analysis/data/derived/seroprevalence/Aden/aden_excess_complete_sero.RDS")


  return(list(fit = fit, fit_complete = fit_complete, fit_complete_sero = fit_complete_sero))

}

# fit with 100000 mcmc iterations
# we suggest running on a HPC as is computationally intensive

aden_excess <- aden_fit_excess_smooth(100000,1)

# vary IFR in model fit























