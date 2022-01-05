
#' Fit squire model using rw splines
#'
#' @inheritParams squire::pmcmc
#' @inheritParams squire::parameters_explicit_SEEIR
#' @param hosp_beds General Hospital Beds
#' @param icu_beds ICU Beds
#' @param rw_duration Random Walk/Spline Duration. Default = 14 days
#' @param initial_pars List of starting mcmc pars. Default = NULL
#'
#' @return Model fit from [squire:::pmcmc]
#' @export
#'
fit_spline_rt <- function(data,
                          country,
                          population,
                          reporting_fraction,
                          n_mcmc = 10000,
                          replicates = 100,
                          rw_duration = 14,
                          hosp_beds = 10000000000,
                          icu_beds = 10000000000,
                          log_likelihood = NULL,
                          initial_pars = NULL) {

  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------

  # order data
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC SETUP
  ## -----------------------------------------------------------------------------

  # dat_0 is just the current date now
  date_0 <- max(data$date)

  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$date[which(data$deaths>0)][1]

  # We set the R0_change here to be 1 everywhere to effectively turn off mobility
  R0_change <- rep(1, nrow(data))
  date_R0_change <- data$date
  R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
  date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]

  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  n_chains <- 3 # number of chains
  start_adaptation <- as.integer(0.1*n_mcmc) # how long before adapting

  # parallel call
  suppressWarnings(future::plan(future::multiprocess()))

  # Defualt parameter edges for pmcmc
  R0_min <- 1.6
  R0_max <- 5.6
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-30

  # These 4 parameters do nothign as setting R0_change to 1
  Meff_min <- -2
  Meff_max <- 2
  Meff_pl_min <- 0
  Meff_pl_max <- 1
  Rt_shift_min <- 0
  Rt_shift_max <- 0.001
  Rt_shift_scale_min <- 0.1
  Rt_shift_scale_max <- 10


  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------

  date_start <- data$date[which(cumsum(data$deaths)>10)[1]] - 30
  R0_start <- 3

  # These are the the initial conditions now loaded from our previous run.
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))

  # again these all do nothing
  Meff_start <- min(max(0, Meff_min), Meff_max)
  Meff_pl_start <- min(max(0.5, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(0.0005, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(5, Rt_shift_scale_min), Rt_shift_scale_max)

  # Our random walk parameters start after the Meff change
  # Basically just set this suitably far back in the past
  date_Meff_change <- date_start - 1

  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  last_shift_date <- as.Date(date_Meff_change) + 1
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place

  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks.
  rw_needed <- as.numeric(ceiling(remaining_days/Rt_rw_duration))

  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter set up
  ## -----------------------------------------------------------------------------

  # PMCMC Parameters
  pars_init = list('start_date' = date_start,
                   'R0' = R0_start,
                   'Meff' = Meff_start,
                   'Meff_pl' = Meff_pl_start,
                   "Rt_shift" = 0,
                   "Rt_shift_scale" = Rt_shift_scale_start)
  pars_min = list('start_date' = first_start_date,
                  'R0' = R0_min,
                  'Meff' = Meff_min,
                  'Meff_pl' = Meff_pl_min,
                  "Rt_shift" = Rt_shift_min,
                  "Rt_shift_scale" = Rt_shift_scale_min)
  pars_max = list('start_date' = last_start_date,
                  'R0' = R0_max,
                  'Meff' = Meff_max,
                  'Meff_pl' = Meff_pl_max,
                  "Rt_shift" = Rt_shift_max,
                  "Rt_shift_scale" = Rt_shift_scale_max)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE,
                       'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

  # add in the spline list
  pars_init <- append(pars_init, pars_init_rw)
  pars_min <- append(pars_min, pars_min_rw)
  pars_max <- append(pars_max, pars_max_rw)
  pars_discrete <- append(pars_discrete, pars_discrete_rw)

  # use initial pars if provided
  if (!is.null(initial_pars)) {
    pars_init <- initial_pars
  }

  # Covriance Matrix
  proposal_kernel <- diag(length(names(pars_init))) * 0.3
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
  proposal_kernel["start_date", "start_date"] <- 1.5

  # MCMC Functions - Prior and Likelihood Calculation
  logprior <- function(pars){
    ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
      dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
      dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE)

    # get rw spline parameters
    if(any(grepl("Rt_rw", names(pars)))) {
      Rt_rws <- pars[grepl("Rt_rw", names(pars))]
      for (i in seq_along(Rt_rws)) {
        ret <- ret + dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
      }
    }
    return(ret)
  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------

  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)

  # run the pmcmc
  res <- squire::pmcmc(data = data,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = log_likelihood,
                       reporting_fraction = reporting_fraction,
                       squire_model = squire:::deterministic_model(),
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       proposal_kernel = proposal_kernel,
                       population = population,
                       baseline_contact_matrix = mix_mat,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       Rt_args = squire:::Rt_args_list(
                         date_Meff_change = date_Meff_change,
                         scale_Meff_pl = TRUE,
                         Rt_shift_duration = 1,
                         Rt_rw_duration = Rt_rw_duration),
                       burnin = ceiling(n_mcmc/10),
                       seeding_cases = 5,
                       replicates = replicates,
                       required_acceptance_ratio = 0.20,
                       start_adaptation = start_adaptation,
                       baseline_hosp_bed_capacity = hosp_beds,
                       baseline_ICU_bed_capacity = icu_beds)


  ## remove things so they don't atke up so much memory when you save them :)

  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  # remove states to keep object memory save down
  for(i in seq_along(res$pmcmc_results$chains)) {
    res$pmcmc_results$chains[[i]]$states <- NULL
    res$pmcmc_results$chains[[i]]$covariance_matrix <- tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
  }

  return(res)

}

#' Fit squire model using rw splines
#'
#' @inheritParams squire::pmcmc
#' @inheritParams squire::parameters_explicit_SEEIR
#' @inheritParams fit_spline_rt
#' @param prob_death_scale Probability of death given infection scalar.
#'
#' @return Model fit from [squire:::pmcmc]
#' @export
#'
fit_spline_rt_ifr_change <- function(data,
                                     country,
                                     population,
                                     reporting_fraction,
                                     log_likelihood = NULL,
                                     n_mcmc = 10000,
                                     replicates = 100,
                                     rw_duration = 14,
                                     hosp_beds = 10000000000,
                                     icu_beds = 10000000000,
                                     prob_death_scale) {

  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------

  # order data
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC SETUP
  ## -----------------------------------------------------------------------------

  # dat_0 is just the current date now
  date_0 <- max(data$date)

  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$date[which(data$deaths>0)][1]

  # We set the R0_change here to be 1 everywhere to effectively turn off mobility
  R0_change <- rep(1, nrow(data))
  date_R0_change <- data$date
  R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
  date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]

  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  n_chains <- 3 # number of chains
  start_adaptation <- as.integer(0.1*n_mcmc) # how long before adapting

  # parallel call
  suppressWarnings(future::plan(future::multiprocess()))

  # Defualt parameter edges for pmcmc
  R0_min <- 1.6
  R0_max <- 5.6
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-30

  # These 4 parameters do nothign as setting R0_change to 1
  Meff_min <- -2
  Meff_max <- 2
  Meff_pl_min <- 0
  Meff_pl_max <- 1
  Rt_shift_min <- 0
  Rt_shift_max <- 0.001
  Rt_shift_scale_min <- 0.1
  Rt_shift_scale_max <- 10


  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------

  date_start <- data$date[which(cumsum(data$deaths)>10)[1]] - 30
  R0_start <- 3

  # These are the the initial conditions now loaded from our previous run.
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))

  # again these all do nothing
  Meff_start <- min(max(0, Meff_min), Meff_max)
  Meff_pl_start <- min(max(0.5, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(0.0005, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(5, Rt_shift_scale_min), Rt_shift_scale_max)

  # Our random walk parameters start after the Meff change
  # Basically just set this suitably far back in the past
  date_Meff_change <- date_start - 1

  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  last_shift_date <- as.Date(date_Meff_change) + 1
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place

  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks.
  rw_needed <- as.numeric(ceiling(remaining_days/Rt_rw_duration))

  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter set up
  ## -----------------------------------------------------------------------------

  # PMCMC Parameters
  pars_init = list('start_date' = date_start,
                   'R0' = R0_start,
                   'Meff' = Meff_start,
                   'Meff_pl' = Meff_pl_start,
                   "Rt_shift" = 0,
                   "Rt_shift_scale" = Rt_shift_scale_start)
  pars_min = list('start_date' = first_start_date,
                  'R0' = R0_min,
                  'Meff' = Meff_min,
                  'Meff_pl' = Meff_pl_min,
                  "Rt_shift" = Rt_shift_min,
                  "Rt_shift_scale" = Rt_shift_scale_min)
  pars_max = list('start_date' = last_start_date,
                  'R0' = R0_max,
                  'Meff' = Meff_max,
                  'Meff_pl' = Meff_pl_max,
                  "Rt_shift" = Rt_shift_max,
                  "Rt_shift_scale" = Rt_shift_scale_max)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE,
                       'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6)

  # add in the spline list
  pars_init <- append(pars_init, pars_init_rw)
  pars_min <- append(pars_min, pars_min_rw)
  pars_max <- append(pars_max, pars_max_rw)
  pars_discrete <- append(pars_discrete, pars_discrete_rw)

  # use initial pars if provided
  if (!is.null(initial_pars)) {
    pars_init <- initial_pars
  }

  # Covriance Matrix
  proposal_kernel <- diag(length(names(pars_init))) * 0.3
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
  proposal_kernel["start_date", "start_date"] <- 1.5

  # MCMC Functions - Prior and Likelihood Calculation
  logprior <- function(pars){
    ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
      dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
      dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE)

    # get rw spline parameters
    if(any(grepl("Rt_rw", names(pars)))) {
      Rt_rws <- pars[grepl("Rt_rw", names(pars))]
      for (i in seq_along(Rt_rws)) {
        ret <- ret + dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
      }
    }
    return(ret)
  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------

  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)

  # run the pmcmc
  res <- squire::pmcmc(data = data,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = log_likelihood,
                       reporting_fraction = reporting_fraction,
                       squire_model = squire:::deterministic_model(),
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       proposal_kernel = proposal_kernel,
                       population = population,
                       baseline_contact_matrix = mix_mat,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       Rt_args = squire:::Rt_args_list(
                         date_Meff_change = date_Meff_change,
                         scale_Meff_pl = TRUE,
                         Rt_shift_duration = 1,
                         Rt_rw_duration = Rt_rw_duration),
                       burnin = ceiling(n_mcmc/10),
                       seeding_cases = 5,
                       replicates = replicates,
                       required_acceptance_ratio = 0.20,
                       start_adaptation = start_adaptation,
                       baseline_hosp_bed_capacity = hosp_beds,
                       baseline_ICU_bed_capacity = icu_beds,
                       prob_severe_death_treatment = vapply(squire:::default_probs()$prob_severe_death_treatment*prob_death_scale,
                                                            min, numeric(1), 0.95),
                       prob_non_severe_death_treatment = vapply(squire:::default_probs()$prob_non_severe_death_treatment*prob_death_scale,
                                                                min, numeric(1), 0.95))


  ## remove things so they don't atke up so much memory when you save them :)

  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  # remove states to keep object memory save down
  for(i in seq_along(res$pmcmc_results$chains)) {
    res$pmcmc_results$chains[[i]]$states <- NULL
    res$pmcmc_results$chains[[i]]$covariance_matrix <- tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
  }

  return(res)

}

#' @noRd
lusaka_log_likelihood <- function(pars, data, squire_model, model_params, pars_obs, n_particles,
                                forecast_days = 0, return = "ll", Rt_args, interventions, ...) {

  switch(return, full = {
    save_particles <- TRUE
    full_output <- TRUE
    pf_return <- "sample"
  }, ll = {
    save_particles <- FALSE
    forecast_days <- 0
    full_output <- FALSE
    pf_return <- "single"
  }, {
    stop("Unknown return type to calc_loglikelihood")
  })

  squire:::assert_in(c("R0", "start_date"), names(pars), message = "Must specify R0, start date to infer")
  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)
  R0_change <- interventions$R0_change
  date_R0_change <- interventions$date_R0_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change

  if (is.null(date_R0_change)) {
    tt_beta <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                    change = R0_change, start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }

  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(interventions$contact_matrix_set)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change,
                                                         , ]
  }

  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = interventions$ICU_bed_capacity[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }

  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = interventions$hosp_bed_capacity[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }

  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = interventions$max_vaccine[-1], start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt), starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }

  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_vaccine_efficacy_infection <- tt_list$tt
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,
                                                                                       , ]
  }

  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  }
  else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date, steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_vaccine_efficacy_disease <- tt_list$tt
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,
                                                     , ]
  }

  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change, R0 = R0, date_R0_change = date_R0_change,
                                   pars = pars, Rt_args = Rt_args)
  beta_set <- squire:::beta_est(squire_model = squire_model, model_params = model_params,
                                R0 = R0)
  model_params$beta_set <- beta_set

  if (inherits(squire_model, "stochastic")) {
    pf_result <- squire:::run_particle_filter(data = data, squire_model = squire_model,
                                              model_params = model_params, model_start_date = start_date,
                                              obs_params = pars_obs, n_particles = n_particles,
                                              forecast_days = forecast_days, save_particles = save_particles,
                                              full_output = full_output, return = pf_return)
  }
  else if (inherits(squire_model, "deterministic")) {
    pf_result <- run_deterministic_comparison_lusaka(data = data,
                                                   squire_model = squire_model, model_params = model_params,
                                                   model_start_date = start_date, obs_params = pars_obs,
                                                   forecast_days = forecast_days, save_history = save_particles,
                                                   return = pf_return)
  }
  pf_result

}

#' @noRd
run_deterministic_comparison_lusaka <- function(data, squire_model, model_params, model_start_date = "2020-02-02",
                                              obs_params = list(
                                                phi_cases = 0.1,
                                                k_cases = 2,
                                                phi_death = 1,
                                                k_death = 2,
                                                exp_noise = 1e+06
                                              ), forecast_days = 0, save_history = FALSE,
                                              return = "ll") {

  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$date[data$deaths > 0][1], "%Y-%m-%d") <
      as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }

  # set up as normal
  data <- squire:::particle_filter_data(data = data, start_date = model_start_date,
                                        steps_per_day = round(1/model_params$dt))
  # correct for weekly deaths
  data$day_end[nrow(data)] <- data$day_start[nrow(data)] + 7
  data$step_end[nrow(data)] <- data$step_start[nrow(data)] + 7

  # back to normal
  model_params$tt_beta <- round(model_params$tt_beta * model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix *
                                            model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds *
                                       model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds *
                                      model_params$dt)

  # steps as normal
  steps <- c(0, data$day_end)
  fore_steps <- seq(data$day_end[nrow(data)], length.out = forecast_days + 1L)
  steps <- unique(c(steps, fore_steps))

  # run model
  model_func <- squire_model$odin_model(user = model_params,
                                        unused_user_action = "ignore")
  out <- model_func$run(t = seq(0, tail(steps, 1), 1), atol = 1e-6, rtol = 1e-6)
  index <- squire:::odin_index(model_func)

  # get deaths for comparison
  Ds <- diff(rowSums(out[c(data$day_end[2]-7, data$day_end[-1]), index$D]))
  Ds[Ds < 0] <- 0
  deaths <- data$deaths[-1]
  ll <- squire:::ll_nbinom(deaths, Ds, obs_params$phi_death, 1,
                  obs_params$exp_noise)


  # and wrap up as normal
  date <- data$date[[1]] + seq_len(nrow(out)) - 1L
  rownames(out) <- as.character(date)
  attr(out, "date") <- date
  pf_results <- list()
  pf_results$log_likelihood <- sum(ll)
  if (save_history) {
    pf_results$states <- out
  }
  else if (return == "single") {
    pf_results$sample_state <- out[nrow(out), ]
  }
  if (return == "ll") {
    ret <- pf_results$log_likelihood
  }
  else if (return == "sample") {
    ret <- pf_results$states
  }
  else if (return == "single" || return == "full") {
    ret <- pf_results
  }
  ret
}



#' Fit squire model using rw splines
#'
#' @inheritParams squire::pmcmc
#' @inheritParams squire::parameters_explicit_SEEIR
#' @param hosp_beds General Hospital Beds
#' @param icu_beds ICU Beds
#' @param rw_duration Random Walk/Spline Duration. Default = 14 days
#' @param reporting_fraction_bounds Region to vary reporting fraction
#'
#' @return Model fit from [squire:::pmcmc]
#' @export
#'
fit_spline_rf_fit <- function(data,
                          country,
                          population,
                          reporting_fraction=1,
                          reporting_fraction_bounds=NULL,
                          log_likelihood = NULL,
                          n_mcmc = 10000,
                          replicates = 100,
                          rw_duration = 14,
                          hosp_beds = 10000000000,
                          icu_beds = 10000000000,
                          initial_pars = NULL,
                          sero_det = NULL,
                          sero_df = NULL,
                          pcr_det = NULL,
                          pcr_df = NULL,
                          ...) {

  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------

  # order data
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC SETUP
  ## -----------------------------------------------------------------------------

  # dat_0 is just the current date now
  date_0 <- max(data$date)

  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$date[which(data$deaths>0)][1]

  # We set the R0_change here to be 1 everywhere to effectively turn off mobility
  R0_change <- rep(1, nrow(data))
  date_R0_change <- data$date
  R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
  date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]

  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  n_chains <- 3 # number of chains
  start_adaptation <- as.integer(0.1*n_mcmc) # how long before adapting

  # parallel call
  suppressWarnings(future::plan(future::multiprocess()))

  # Defualt parameter edges for pmcmc
  R0_min <- 1.6
  R0_max <- 5.6
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-30

  # These 4 parameters do nothign as setting R0_change to 1
  Meff_min <- -2
  Meff_max <- 2
  Meff_pl_min <- 0
  Meff_pl_max <- 1
  Rt_shift_min <- 0
  Rt_shift_max <- 0.001
  Rt_shift_scale_min <- 0.1
  Rt_shift_scale_max <- 10


  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------

  date_start <- data$date[which(cumsum(data$deaths)>10)[1]] - 30
  R0_start <- 3

  # These are the the initial conditions now loaded from our previous run.
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))

  # again these all do nothing
  Meff_start <- min(max(0, Meff_min), Meff_max)
  Meff_pl_start <- min(max(0.5, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(0.0005, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(5, Rt_shift_scale_min), Rt_shift_scale_max)

  # Our random walk parameters start after the Meff change
  # Basically just set this suitably far back in the past
  date_Meff_change <- date_start - 1

  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  last_shift_date <- as.Date(date_Meff_change) + 1
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place

  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks.
  rw_needed <- as.numeric(ceiling(remaining_days/Rt_rw_duration))

  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter set up
  ## -----------------------------------------------------------------------------

  # PMCMC Parameters
  pars_init = list('start_date' = date_start,
                   'R0' = R0_start,
                   'Meff' = Meff_start,
                   'Meff_pl' = Meff_pl_start,
                   "Rt_shift" = 0,
                   "Rt_shift_scale" = Rt_shift_scale_start)
  pars_min = list('start_date' = first_start_date,
                  'R0' = R0_min,
                  'Meff' = Meff_min,
                  'Meff_pl' = Meff_pl_min,
                  "Rt_shift" = Rt_shift_min,
                  "Rt_shift_scale" = Rt_shift_scale_min)
  pars_max = list('start_date' = last_start_date,
                  'R0' = R0_max,
                  'Meff' = Meff_max,
                  'Meff_pl' = Meff_pl_max,
                  "Rt_shift" = Rt_shift_max,
                  "Rt_shift_scale" = Rt_shift_scale_max)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE,
                       'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = 2, exp_noise = 1e6,
                  sero_det = sero_det, pcr_det = pcr_det, sero_df = sero_df, pcr_df = pcr_df)

  null_obs <- c()
  for(i in seq_along(pars_obs)) {
    null_obs[i] <- is.null(pars_obs[[i]])
  }
  pars_obs <- pars_obs[-which(null_obs)]

  # add in the spline list
  pars_init <- append(pars_init, pars_init_rw)
  pars_min <- append(pars_min, pars_min_rw)
  pars_max <- append(pars_max, pars_max_rw)
  pars_discrete <- append(pars_discrete, pars_discrete_rw)

  # add reporting bounds if given
  if(!is.null(reporting_fraction_bounds)){
    pars_init <- append(pars_init, c("rf"=reporting_fraction_bounds[1]))
    pars_min <- append(pars_min, c("rf"=reporting_fraction_bounds[2]))
    pars_max <- append(pars_max, c("rf"=reporting_fraction_bounds[3]))
    pars_discrete <- append(pars_discrete, c("rf"=FALSE))
  }

  # Covariance Matrix
  proposal_kernel <- diag(length(names(pars_init))) * 0.3
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
  proposal_kernel["start_date", "start_date"] <- 1.5

  # MCMC Functions - Prior and Likelihood Calculation
  logprior <- function(pars){
    ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
      dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
      dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE)

    # get rw spline parameters
    if(any(grepl("Rt_rw", names(pars)))) {
      Rt_rws <- pars[grepl("Rt_rw", names(pars))]
      for (i in seq_along(Rt_rws)) {
        ret <- ret + dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
      }
    }
    return(ret)
  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------

  # mixing matrix - assume is same as country as whole
  mix_mat <- squire::get_mixing_matrix(country)

  # run the pmcmc
  res <- squire::pmcmc(data = data,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = NULL,
                       reporting_fraction = reporting_fraction,
                       squire_model = squire:::deterministic_model(),
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       population = population,
                       baseline_contact_matrix = mix_mat,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       Rt_args = squire:::Rt_args_list(
                         date_Meff_change = date_Meff_change,
                         scale_Meff_pl = TRUE,
                         Rt_shift_duration = 1,
                         Rt_rw_duration = Rt_rw_duration),
                       burnin = ceiling(n_mcmc/10),
                       seeding_cases = 5,
                       replicates = replicates,
                       required_acceptance_ratio = 0.20,
                       start_adaptation = start_adaptation,
                       baseline_hosp_bed_capacity = hosp_beds,
                       baseline_ICU_bed_capacity = icu_beds,
                       ...)


  ## remove things so they don't atke up so much memory when you save them :)

  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  # remove states to keep object memory save down
  for(i in seq_along(res$pmcmc_results$chains)) {
    res$pmcmc_results$chains[[i]]$states <- NULL
    res$pmcmc_results$chains[[i]]$covariance_matrix <- tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
  }

  # Add sero_det and pcr_det to output
  res$pmcmc_results$inputs$pars_obs$sero_det <- sero_det
  res$pmcmc_results$inputs$pars_obs$pcr_det <- pcr_det

  return(res)

}




