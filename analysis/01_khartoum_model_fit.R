# fit model to COVID-19 and scan across reporting fraction ot identify best fit to symptomatic prevalence
library(lubridate)
library(tidyverse)
library(odin)
library(squire)
library(dde)
devtools::load_all()

### fit to lusaka COVID-19 deaths
khartoum_fit <- function(n_mcmc,
                         rf = 1,
                         city_age = "same",
                       hosp_beds = 4640,
                       icu_beds = 75,
                       pop_size = 8877147,
                       initial_pars = NULL){


  # -----------------------------------------------------------------------------
  # Step 1: Country specific demography and mix mat
  # -----------------------------------------------------------------------------

  # Sudan demgraphy
  dem <- squire::get_population("Sudan")

  # city_age <- "older"
  if (city_age == "younger") {
    altered_pop <- dem$n * (dexp(1:17, 1/50)/mean(dexp(1:17, 1/500)))
    altered_pop <- altered_pop / sum(altered_pop)
    population <- round(altered_pop*pop_size)
  } else if (city_age == "older") {
    altered_pop <- dem$n * rev((dexp(1:17, 1/50))/mean(dexp(1:17, 1/50)))
    altered_pop <- altered_pop / sum(altered_pop)
    population <- round(altered_pop*pop_size)
  } else {
    population <- round((dem$n/sum(dem$n))*pop_size)
  }

  # matrix
  mix_mat <- squire::get_mixing_matrix("Sudan")

  # -----------------------------------------------------------------------------
  # Step 2: Soruces for other key parameters
  # -----------------------------------------------------------------------------

  # ------------------------------
  # BEDS
  # ------------------------------

  # icu bed sources

  # mid April/May Covid treatament in icus were provided in a few gvt facilities. This included icu access.
  # Anecdotally though there is a sense that this was at peak in April May June and since have started to come down
  # https://go.ifrc.org/reports/13075 # gives a live upper estimate
  # https://go.ifrc.org/reports/13042
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132516/ - "74-bed ICUs are almost fully occupied year-round"

  # Hospital bed sources

  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6132516/ - 1737 across the 4 government tertiary care hospitals
  # per capita reference https://sudannextgen.com/wp-content/uploads/2019/10/Statistics-on-Hospital-sector-in-Sudan-2.pdf has much higher
  # 4640

  # get our deaths source
  deaths <- readRDS("analysis/data/derived/deaths_time_series/khartoum_covid_deaths.rds")

  # suitable starting pars
  if(is.null(initial_pars)) {
    initial_pars <- readRDS("analysis/data/derived/model_fits/Lusaka/pars_init.rds")
  }

  fit <- fit_spline_rt_khartoum(data=deaths,
                       country="Sudan",
                       population=population,
                       reporting_fraction=rf,
                       n_mcmc = n_mcmc,
                       replicates = 10,
                       rw_duration = 14,
                       hosp_beds = hosp_beds,
                       icu_beds = icu_beds,
                       initial_pars = initial_pars)

  saveRDS(fit, paste0("analysis/data/derived/model_fits/Khartoum/khartoum_fit_rf_", rf, ".rds"))

  return(fit)

}



# we suggest running on a HPC as is computationally intensive
n_mcmc <- 50000
pars_init <- readRDS("analysis/data/derived/model_fits/Khartoum/pars_init.rds")
khartoum_model_fits <- list()

rfs <- c(seq(0.02,0.06, 0.01), 0.07, 0.08, 0.09, 0.1,0.2,0.4, 0.045)
for(i in seq_along(rfs)) {
  khartoum_model_fits[[i]] <- khartoum_fit(n_mcmc, rfs[i], initial_pars = pars_init[[i]])
}

# function to calculate the likelihood of the fits against symptomatic prevalence source
khartoum_ll <- function(res){

  date_0 <- max(res$pmcmc_results$inputs$data$date)

  ## Symptomatic infection prevalence estimate as of 2nd June which we compare to cum_hosp_inc + cum_icu_inc
  # https://www.researchgate.net/publication/342122988_Estimation_of_Coronavirus_COVID-19_Infections_in_Khartoum_State_Sudan
  p_e_low <- 0.083
  p_e_med <- 0.11
  p_e_high <- 0.137

  # first lets get symptomatic totals for individuals sufficiently severe to need treatment beofre 2nd June and older than 15
  hosp_inc <- squire::format_output(res, c("hospital_incidence"), date_0 = date_0, reduce_age = FALSE) %>%
    filter(age_group > 3)

  icu_inc <- squire::format_output(res, c("ICU_incidence"), date_0 = date_0, reduce_age = FALSE) %>%
    filter(age_group > 3)

  # and then also total infs
  infs <- squire::format_output(res, c("infections"), date_0 = date_0, reduce_age = FALSE) %>%
    filter(age_group > 3)

  # symptomatic ratios
  symptoms <- read.csv("analysis/data/raw/symptom_probability.csv")
  age_matches <- unlist(lapply(sapply(sort(unique(infs$age_group)), grep, symptoms$age_group), "[[", 1))
  ages_matched <- sort(unique(infs$age_group))

  # we will use these probabilities of symptoms in all cases to calculate probability
  # of symptoms in non hospitalised patients:
  symp_med <- symptoms$mid - res$parameters$prob_hosp
  symp_low <- symptoms$low - res$parameters$prob_hosp
  symp_high <- symptoms$high - res$parameters$prob_hosp

  # estimate of symptomatic numbers per age
  infs$non_severe <- infs$y - hosp_inc$y - icu_inc$y
  infs$symptomatic <- hosp_inc$y + icu_inc$y + infs$non_severe*symp_med[age_matches[match(infs$age_group, ages_matched)]]
  infs$symptomatic_low <- hosp_inc$y + icu_inc$y + infs$non_severe*symp_low[age_matches[match(infs$age_group, ages_matched)]]
  infs$symptomatic_high <- hosp_inc$y + icu_inc$y + infs$non_severe*symp_high[age_matches[match(infs$age_group, ages_matched)]]

  # summarise back up to date
  infs_sum <- infs %>% group_by(replicate, t, date) %>%
    summarise(symptomatic = sum(symptomatic, na.rm = TRUE),
              symptomatic_low = sum(symptomatic_low, na.rm = TRUE),
              symptomatic_high = sum(symptomatic_high, na.rm = TRUE))
  infs_sum$p_e <- infs_sum$symptomatic / sum(res$parameters$population[-c(1:3)])
  infs_sum$p_e_low <- infs_sum$symptomatic_low / sum(res$parameters$population[-c(1:3)])
  infs_sum$p_e_high <- infs_sum$symptomatic_high / sum(res$parameters$population[-c(1:3)])

  infs_sum <- infs_sum %>% group_by(replicate) %>%
    mutate(cumu_symp_pe_med = cumsum(p_e),
           cumu_symp_pe_low = cumsum(p_e_low),
           cumu_symp_pe_high = cumsum(p_e_high))

  # how many cases would this be
  expected_cases_low <- round(p_e_low * sum(res$parameters$population))
  expected_cases_med <- round(p_e_med * sum(res$parameters$population))
  expected_cases_high <- round(p_e_high * sum(res$parameters$population))

  # how many were model observed
  observed_cases_low <- round(infs_sum$cumu_symp_pe_low[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))
  observed_cases_med <- round(infs_sum$cumu_symp_pe_med[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))
  observed_cases_high <- round(infs_sum$cumu_symp_pe_high[which(infs_sum$date == "2020-06-02")] * sum(res$parameters$population))

  # get likelihoods against cumulative cases described by a poisson

  ll_cases_low <- dpois(observed_cases_low, lambda = expected_cases_low, log = TRUE)
  ll_cases_med <- dpois(observed_cases_med, lambda = expected_cases_med, log = TRUE)
  ll_cases_high <- dpois(observed_cases_high, lambda = expected_cases_high, log = TRUE)

  data.frame("ll_low" = ll_cases_low,
                   "ll_med" = ll_cases_med,
                   "ll_high" = ll_cases_high,
             "p_e_low" = observed_cases_low/sum(res$parameters$population),
             "p_e_med" = observed_cases_med/sum(res$parameters$population),
             "p_e_high" = observed_cases_high/sum(res$parameters$population),
                   "replicate" = seq_along(ll_cases_low),
                   "rf" = res$pmcmc_results$inputs$pars_obs$phi_death
  )

}

ll <- do.call(rbind, lapply(khartoum_model_fits, khartoum_ll))
saveRDS(ll, "analysis/data/derived/seroprevalence/Khartoum/likelihood_rf.rds")
