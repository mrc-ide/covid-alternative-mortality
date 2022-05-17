# fit model to COVID-19 and excess mortality data from Addis Ababa and estimate seroprevalence

library(lubridate)
library(tidyverse)
library(odin)
library(squire)
library(dde)
devtools::load_all()

### fit to lusaka COVID-19 deaths

lusaka_downstream <- function(fit, rf, death_source){

  # from model fits to Van Elslande May/June Fits
  sero_det <- sero_det("igg", igg_sens = 0.95, igg_conv = 13.3)
  sero_fit <- seroprev_df_det(fit, sero_det = sero_det)

  # from Hay et al 2021 Science (actually from preprint)
  pcr_sens <- 0.95
  pcr_det <- c(9.206156e-13, 9.206156e-13, 3.678794e-01, 9.645600e-01,
               9.575796e-01, 9.492607e-01, 9.393628e-01, 9.276090e-01,
               9.136834e-01, 8.972309e-01, 8.778578e-01, 8.551374e-01,
               8.286197e-01, 7.978491e-01, 7.623916e-01, 7.218741e-01,
               6.760375e-01, 6.248060e-01, 5.683688e-01, 5.072699e-01,
               4.525317e-01, 4.036538e-01, 3.600134e-01, 3.210533e-01,
               2.862752e-01, 2.552337e-01, 2.275302e-01, 2.028085e-01,
               1.807502e-01, 1.610705e-01, 1.435151e-01, 1.278563e-01,
               1.138910e-01, 1.014375e-01, 9.033344e-02)
  pcr_det <- (pcr_det/max(pcr_det))*pcr_sens
  pcr_det <- c(pcr_det, rep(0, length(sero_det)-length(pcr_det)))

  # combined pcr and sero detectability from infection
  # so we lag the sero_detectability (defined from symptom onset) by 4 days
  combined_det <-   1-((1-c(0, 0, 0, 0, head(sero_det, -4))) * (1-pcr_det))

  sero_fit <- sero_fit %>%
    group_by(replicate) %>%
    mutate(pcr_positive = roll_func(.data$infections, pcr_det),
           pcr_perc = .data$pcr_positive/max(.data$S,na.rm = TRUE),
           combined_positive = roll_func(.data$infections, combined_det),
           combined_perc = .data$combined_positive/max(.data$S,na.rm = TRUE)) %>%
    ungroup

  ### Calculate Log-likelihood
  log_lik_sero_fun <- function(sero_fit){

    sero_fit <- na.omit(sero_fit)

    N_lus <- 2687
    sero_pred_1 <- sero_fit$sero_perc[sero_fit$date==as.Date("2020-07-11")]
    ll_sero_15 <- dbinom(x = as.integer((2.1/10.6 * 9.1)/100*N_lus), size = N_lus, prob = sero_pred_1, log = TRUE)

    sero_pred_2 <- sero_fit$sero_perc[sero_fit$date==as.Date("2020-07-25")]
    ll_sero_29 <- dbinom(x = as.integer(((2.1/10.6 * 9.1) + (7.6/10.6 * 9.1))/100*N_lus), size = N_lus, prob = sero_pred_2, log = TRUE)

    pcr_pred_1 <- sero_fit$pcr_perc[sero_fit$date==as.Date("2020-07-11")]
    ll_pcr_15 <- dbinom(x = as.integer((7.6/10.6 * 9.1)/100*N_lus), size = N_lus, prob = pcr_pred_1, log = TRUE)

    combined_pred <- sero_fit$combined_perc[sero_fit$date==as.Date("2020-07-11")]
    ll_combined_15 <- dbinom(x = as.integer(0.091*N_lus), size = N_lus, prob = combined_pred, log = TRUE)

    return(data.frame("replicate" = seq_along(ll_sero_15),
                      "sero_15" = ll_sero_15,
                      "sero_29" = ll_sero_29,
                      "pcr_15" = ll_pcr_15,
                      "combined_15" = ll_combined_15))
  }

  ll <- log_lik_sero_fun(sero_fit)

  # add these here for ease
  ll$death_source <- death_source
  ll$rf <- rf
  sero_fit$death_source <- death_source
  sero_fit$rf <- rf

  res <- list(fit = fit, sero = sero_fit, ll = ll)
  return(res)
}

lusaka_fit <- function(n_mcmc, rf = 1, death_source = "total_burial_registry_comp",
                       hosp_beds = 1e10, icu_beds = 1e10, lusaka_pop_size = 3360183){

  if(!death_source %in% c("total_burial_registry_comp",
                          "total_burial_registry_comp_mid",
                          "total_burial_registry_comp_strict",
                          "total_sample_effort",
                          "total_sample_effort_mid",
                          "total_sample_effort_strict",
                          "official")) {
    stop("Death source not found in data")
  }


  ## Population size for Lusaka Province
  # 2010: 2191225, 2222812 or 2238569 (https://en.wikipedia.org/wiki/Lusaka_Province)
  # 2014: 2669249 (https://zambia.opendataforafrica.org/apps/atlas/Lusaka)
  # 2018: 3186336 (https://en.wikipedia.org/wiki/Provinces_of_Zambia)
  # 2019: 3308438 (https://www.citypopulation.de/en/zambia/admin/05__lusaka/)
  # 2020: 3360183 (https://zambia.opendataforafrica.org/thrqjfb/population-and-demographic-projections-2011-2035?regionId=ZM-09)
  # 2021: 3484394 (https://zambia.opendataforafrica.org/thrqjfb/population-and-demographic-projections-2011-2035?regionId=ZM-09)

  # 2020 Age distribution for Lusaka Province from opendataforafrica scales for 2020 total
  # (https://zambia.opendataforafrica.org/thrqjfb/population-and-demographic-projections-2011-2035?regionId=ZM-09)
  lusaka_pop <- c(549475,448008,379841,339600,317148,305521,271318,232188,172365,125531,75157,51883,35587,22219,14960,8570,10812)
  lusaka_pop <- round((lusaka_pop/sum(lusaka_pop)) * lusaka_pop_size)

  # get our deaths source
  if (death_source == "official") {
    deaths <- readRDS("analysis/data/derived/deaths_time_series/lusaka_official_deaths.rds")
    lusaka_data <- deaths
    log_likelihood <- NULL
    starting_pars <- NULL
  } else {
    deaths <- readRDS("analysis/data/derived/deaths_time_series/lusaka_postmortem_deaths.rds")
    lusaka_data <- deaths[ ,c("date", death_source)]
    names(lusaka_data) <- c("date", "deaths")
    log_likelihood <- lusaka_log_likelihood
    starting_pars <- list(
      "start_date" = as.Date("2020-04-27"),
      "R0" = 3.982342,
      "Meff" = -0.4017083,
      "Meff_pl" = 0.5332879,
      "Rt_shift" = 0.0004747717,
      "Rt_shift_scale" = 9.03734,
      "Rt_rw_1" = 0.3629854,
      "Rt_rw_2" = 0.5454458,
      "Rt_rw_3" = 0.5609502,
      "Rt_rw_4" = 0.08409054,
      "Rt_rw_5" = 0.3667386,
      "Rt_rw_6" = 0.03206019,
      "Rt_rw_7" = -0.01544258,
      "Rt_rw_8" = 0.1349962,
      "Rt_rw_9" = 0.08671713
    )
  }

  # hosp bed https://cddep.org/wp-content/uploads/2020/05/National-estimates-of-critical-care-capacity-in-54-African-countries..pdf
  # hosp_beds <- as.integer(3360183/18383956 * (18383956/1000)*2)
  # icu bed https://path.azureedge.net/media/documents/Biomedical_Equipment_for_COVID-19_Zambia_facility_assessment_final_4-20-2021.pdf
  # icu_beds <- 43
  fit <- fit_spline_rt(data=lusaka_data,
                       country="Zambia",
                       population=lusaka_pop,
                       reporting_fraction=rf,
                       n_mcmc = n_mcmc,
                       replicates = 10,
                       rw_duration = 14,
                       log_likelihood = log_likelihood,
                       hosp_beds = hosp_beds,
                       icu_beds = icu_beds,
                       initial_pars = starting_pars)

  saveRDS(fit, paste0("analysis/data/derived/model_fits/Lusaka/lusaka_", death_source, "_rf_", rf, "_fit.rds"))

  res <- lusaka_downstream(fit, rf, death_source)
  saveRDS(res, paste0("analysis/data/derived/seroprevalence/Lusaka/lusaka_", death_source, "_rf_", rf, "_ll.rds"))

  return(res)

}



# fit with 50000 mcmc iterations
# we suggest running on a HPC as is computationally intensive
n_mcmc <- 50000
death_sources <- c("total_burial_registry_comp",
                   "total_burial_registry_comp_mid",
                  "total_burial_registry_comp_strict",
                  "total_sample_effort",
                  "total_sample_effort_mid",
                  "total_sample_effort_strict")

lusaka_mortem_fits <- list()
for(i in seq_along(death_sources)) {
  lusaka_mortem_fits[[death_sources[[i]]]] <- lusaka_fit(n_mcmc, 1, death_sources[i])
}

lusaka_official_fits <- list()
rfs <- seq(0.1, 0.5, 0.1)
for(i in seq_along(rfs)) {
  lusaka_official_fits[[paste0("rf_", rfs[[i]])]] <- lusaka_fit(n_mcmc, rfs[[i]], "official")
}



# RF FIT

lusaka_rf_fit <- function(n_mcmc, rf = 0.3,
                          hosp_beds = 1e10, icu_beds = 1e10,
                          lusaka_pop_size = 3360183){

  death_source <- "official"
  lusaka_pop <- c(549475,448008,379841,339600,317148,305521,271318,232188,172365,125531,75157,51883,35587,22219,14960,8570,10812)
  lusaka_pop <- round((lusaka_pop/sum(lusaka_pop)) * lusaka_pop_size)

  # get our deaths source
    deaths <- readRDS("analysis/data/derived/deaths_time_series/lusaka_official_deaths.rds")
    lusaka_data <- deaths
    log_likelihood <- NULL
    starting_pars <- list(
      "start_date" = as.Date("2020-04-27"),
      "R0" = 3.982342,
      "Meff" = -0.4017083,
      "Meff_pl" = 0.5332879,
      "Rt_shift" = 0.0004747717,
      "Rt_shift_scale" = 9.03734,
      "Rt_rw_1" = 0.3629854,
      "Rt_rw_2" = 0.5454458,
      "Rt_rw_3" = 0.5609502,
      "Rt_rw_4" = 0.08409054,
      "Rt_rw_5" = 0.3667386,
      "Rt_rw_6" = 0.03206019,
      "Rt_rw_7" = -0.01544258,
      "Rt_rw_8" = 0.1349962,
      "Rt_rw_9" = 0.08671713
    )

    # from model fits to Van Elslande May/June Fits
    sero_det <- sero_det("igg", igg_sens = 0.95, igg_conv = 13.3)

    # from Hay et al 2021 Science (actually from preprint)
    pcr_sens <- 0.95
    pcr_det <- c(9.206156e-13, 9.206156e-13, 3.678794e-01, 9.645600e-01,
                 9.575796e-01, 9.492607e-01, 9.393628e-01, 9.276090e-01,
                 9.136834e-01, 8.972309e-01, 8.778578e-01, 8.551374e-01,
                 8.286197e-01, 7.978491e-01, 7.623916e-01, 7.218741e-01,
                 6.760375e-01, 6.248060e-01, 5.683688e-01, 5.072699e-01,
                 4.525317e-01, 4.036538e-01, 3.600134e-01, 3.210533e-01,
                 2.862752e-01, 2.552337e-01, 2.275302e-01, 2.028085e-01,
                 1.807502e-01, 1.610705e-01, 1.435151e-01, 1.278563e-01,
                 1.138910e-01, 1.014375e-01, 9.033344e-02)
    pcr_det <- (pcr_det/max(pcr_det))*pcr_sens
    pcr_det <- c(pcr_det, rep(0, length(sero_det)-length(pcr_det)))

    # combined pcr and sero detectability from infection
    # so we lag the sero_detectability (defined from symptom onset) by 4 days
    combined_det <-   1-((1-c(0, 0, 0, 0, head(sero_det, -4))) * (1-pcr_det))

  fit <- fit_spline_rf_fit(data=lusaka_data,
                       country="Zambia",
                       population=lusaka_pop,
                       reporting_fraction=rf,
                       reporting_fraction_bounds=c(rf, 0.1, 0.7),
                       n_mcmc = n_mcmc,
                       replicates = 10,
                       rw_duration = 14,
                       log_likelihood = log_likelihood,
                       hosp_beds = 1e10,
                       icu_beds = 1e10,
                       initial_pars = starting_pars,
                       pcr_det = combined_det,
                       pcr_df = data.frame(
                         "date_start" = as.Date(c("2020-07-03")),
                         "date_end"  = as.Date(c("2020-07-19")),
                         "pcr_pos"  = 244,
                         "samples"  = 2687
                       ))

  saveRDS(fit, paste0("analysis/data/derived/model_fits/Lusaka/lusaka_", death_source, "_rfest_fit.rds"))

  res <- lusaka_downstream(fit, rf, death_source)
  saveRDS(res, paste0("analysis/data/derived/seroprevalence/Lusaka/lusaka_", death_source, "_rfest_ll.rds"))

  return(res)

}


lusaka_rf_fit <- lusaka_rf_fit(10000)



saveRDS(lusaka_mortem_fits, "analysis/data/derived/seroprevalence/Lusaka/Lusaka_postmortem_lls.rds")
saveRDS(lusaka_official_fits, "analysis/data/derived/seroprevalence/Lusaka/Lusaka_official_lls.rds")
saveRDS(lusaka_rf_fit, "analysis/data/derived/seroprevalence/Lusaka/Lusaka_official_rf_est_ll.rds")
