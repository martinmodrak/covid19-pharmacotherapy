prior <- list(
  observation_model_prior = array(c(10,2,2,10), c(2,2))
)

simulator <- function(N_patients, N_obs_per_patient, N_time, prior) {
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
    #first_time <- 1
    patient_times <- sort(sample(first_time:(first_time + 2 * N_obs_per_patient - 1), N_obs_per_patient))
    times[patients == p] <- patient_times
    
    treatment_start_times[p] = rdunif(1, a = first_time, b = max(times) - 1)
  }
  
  
  observations <- array(0L, c(N_patients, N_time))

  for(p in 1:N_patients) {  
    state <- s_ill
    for(t in 1:N_time) {
      observations[p, t] <- sample(1:N_obs_types, size = 1, prob = observation_model[state,])
      if(state == s_ill) {
        if(runif(1) < p_ill_to_healthy) {
          state <- s_healthy
        }
      }
    }
  }
  
  o_types = array(NA_integer_, N_obs)
  for(n in 1:N_obs) {
    o_types[n] = observations[patients[n], times[n]]
  }
  
  list(
    observed = c(prior, list(
      N_patients = N_patients,
      N_obs = N_obs,
      patients = patients,
      times = times,
      o_types = o_types
    )),
    true = list(
      p_ill_to_healthy = p_ill_to_healthy,
      sensitivity = sensitivity,
      specificity = specificity
      #observation_model = observation_model
    )
  )
}