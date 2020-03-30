prior <- list(
  observation_model_prior = array(c(10,2,2,10), c(2,2)),
  state_intercept_prior_sd = 2,
  viral_load_intercept_prior_sd = 5,
  viral_load_intercept_prior_mean = 20,
  viral_load_sigma_prior_sd = 5,
  
  fixed_prior_sd = 1,
  time_effect_prior_sd = 0.5
)

simulator <- function(N_patients, N_treatments, N_obs_per_patient, N_time, prior) {
  N_obs_types <- 2
  o_neg <-  1
  o_pos <-  2
  
  N_states <- 2
  s_healthy <- 1
  s_ill <- 2
  
  observation_model <- matrix(0, nrow = N_states, ncol = N_obs_types)
  # for(s in 1:N_states) {
  #   observation_model[s,] <- MCMCpack::rdirichlet(1, prior$observation_model_prior[s,])
  # }
  sensitivity <- runif(1, 0.5, 1)
  specificity <- runif(1, 0.5, 1);
  observation_model[s_healthy, o_neg] <- specificity
  observation_model[s_healthy, o_pos] <- 1 - specificity
  observation_model[s_ill, o_neg] <- 1 - sensitivity
  observation_model[s_ill, o_pos] <- sensitivity
  
  p_ill_to_healthy <- runif(1)
  
  N_obs <- N_patients * N_obs_per_patient
  patients <- rep(1:N_patients, each = N_obs_per_patient)
  
  times <- integer(N_obs)
  treatment_start_times <- integer(N_patients)
  for(p in 1:N_patients) {
    first_time <- rdunif(1, a = 1, b = N_time -  2 * N_obs_per_patient + 1)
    patient_times <- sort(sample(first_time:(first_time + 2 * N_obs_per_patient - 1), N_obs_per_patient))
    times[patients == p] <- patient_times
    
    treatment_start_times[p] = rdunif(1, a = first_time, b = max(times) - 1)
  }
  
  patients_per_treatment <- round(N_patients / (N_treatments + 1))
  patients_no_treatment <- N_patients - patients_per_treatment * N_treatments
  if(N_treatments > 0) {
    N_fixed <- N_treatments
    X <- array(0, c(N_obs, N_fixed))

    fixed_prior_sd <- array(c(prior$intercept_sd), N_fixed)
    
    treatment_per_patient <- integer(N_patients)
    treatment_per_patient[1:patients_no_treatment] <- 0
    for(t in 1:N_treatments) {
      start <- (patients_no_treatment + (t - 1) * patients_per_treatment) + 1
      patients_for_treatment <- start:(start + patients_per_treatment - 1)
      treatment_per_patient[patients_for_treatment] <- t
      observation_indices <- patients %in% patients_for_treatment
      X[observation_indices, t] <- 
        observation_times[observation_indices] - treatment_start_times[observation_patients[observation_indices]]
    }
  } else {
    ## All zeroes
    treatment_per_patient <- integer(N_patients)
  }
  
  beta <- rnorm(N_fixed, 0, sd = fixed_prior_sd)
  time_effect <- rnorm(1, 0, sd = prior$time_effect_prior_sd)
  time_effect_shift <- -N_time / 2
  
  o_types_full <- array(0L, c(N_patients, N_time))

  for(p in 1:N_patients) {  
    state <- s_ill
    for(t in 1:N_time) {
      o_types_full[p, t] <- sample(1:N_obs_types, size = 1, prob = observation_model[state,])
      if(state == s_ill) {
        linpred <- sum(X[p,] * beta) + time_effect * (t + time_effect_shift)
        p_ill_to_healthy <- 1 / (1 + exp(-linpred))
        if(runif(1) < p_ill_to_healthy) {
          state <- s_healthy
        }
      }
    }
  }
  
  o_types = array(NA_integer_, N_obs)
  for(n in 1:N_obs) {
    o_types[n] = o_types_full[patients[n], times[n]]
  }
  
  list(
    observed = c(prior, list(
      N_patients = N_patients,
      N_obs = N_obs,
      patients = patients,
      times = times,
      o_types = o_types,
      N_fixed = N_fixed,
      X = X,
      fixed_prior_sd = fixed_prior_sd,
      time_effect_shift = time_effect_shift
    )),
    true = list(
      beta = beta,
      time_effect = time_effect,
      sensitivity = sensitivity,
      specificity = specificity
      #observation_model = observation_model
    )
  )
}