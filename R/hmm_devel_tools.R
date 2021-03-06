sim_hmm <- function(N_patients, N_patients_unknown_load, N_treatments, N_obs_per_patient, 
                    N_time, N_ill_states, prior, N_unknown_shift = 0, max_observation_shift = 0, use_severe_state = TRUE,
                      time_effect = FALSE) {

  # Setup the HMM part
  o_neg <-  1
  o_pos <-  2
  o_severe <- 3
  
  if(use_severe_state) {
    N_states <- N_ill_states + 2
    s_severe <- N_states
  } else {
    N_states <- N_ill_states + 1
    s_severe <- -1
  }
  s_healthy <- 1
  ill_states_shift <- 1
  
  central_ill_state <- as.integer(ceiling(N_ill_states / 2))
  
  state_viral_load_bounds <- seq(from = 0, to = 30, length.out = N_ill_states + 1)
  ill_mean_viral_load <- 0.5 * (state_viral_load_bounds[1 : N_ill_states] + state_viral_load_bounds[2 : (N_ill_states + 1)])
  dim(ill_mean_viral_load) <- N_ill_states

  sensitivity <- runif(1, 0.5, 1)
  specificity <- runif(1, 0.5, 1);
  # observation_model[s_healthy, o_neg] <- specificity
  # observation_model[s_healthy, o_pos] <- 1 - specificity
  # observation_model[s_ill, o_neg] <- 1 - sensitivity
  # observation_model[s_ill, o_pos] <- sensitivity

  observation_sigma <- abs(rnorm(1, sd = prior$observation_sigma_prior_sd))
  transition_thresholds <- sort(rnorm(N_states - 1, 0, prior$transition_thresholds_prior_sd))
  dim(transition_thresholds) <- N_states - 1
  
  state_intercepts_high <- sort(abs(rnorm(N_ill_states - central_ill_state, 0, sd = prior$state_intercept_prior_sd)))
  dim(state_intercepts_high) <- N_ill_states - central_ill_state
  neg_state_intercepts_low <- sort(abs(rnorm(central_ill_state - 1, 0, sd = prior$state_intercept_prior_sd)))
  dim(neg_state_intercepts_low) <- central_ill_state - 1
  state_intercepts <- c(rev(-neg_state_intercepts_low), 0, state_intercepts_high)

  # Setup shared elements
  if(N_unknown_shift > 0) {
    unknown_shift_patients <- sample(1:N_patients, size = N_unknown_shift)
    unknown_shift_ids <- array(0, N_patients)
    unknown_shift_ids[unknown_shift_patients] <- 1:N_unknown_shift
    unknown_shift <- array(NA_integer_, N_unknown_shift)
  } else {
    unknown_shift_patients <- integer(0)
    unknown_shift_ids <- integer(0)
    unknown_shift <- integer(0)
  }
  
  N_obs <- N_patients * N_obs_per_patient
  patients <- rep(1:N_patients, each = N_obs_per_patient)
  
  times <- integer(N_obs)
  treatment_start_times <- integer(N_patients)
  treatment_start_times_observed <- integer(N_patients)
  for(p in 1:N_patients) {
    if(p %in% unknown_shift_patients) {
      first_time <- rdunif(1, a = 1, b = max_observation_shift - 1)
      unknown_shift[unknown_shift_ids[p]] <- first_time - 1
      max_time <- min(first_time + 2 * N_obs_per_patient - 1, N_time - max_observation_shift)
    } else {
      first_time <- rdunif(1, a = 1, b = N_time -  2 * N_obs_per_patient + 1)
      max_time <- first_time + 2 * N_obs_per_patient - 1
    }
    patient_times <- sort(sample(first_time:max_time, N_obs_per_patient))
    times[patients == p] <- patient_times
    
    treatment_start_times[p] = rdunif(1, a = first_time, b = first_time + 3)
    
    
    if(p %in% unknown_shift_patients) {
      treatment_start_times_observed[p] = treatment_start_times[p] - unknown_shift[unknown_shift_ids[p]]
    } else {
      treatment_start_times_observed[p] = treatment_start_times[p]
    }
  }
  
  patients_per_treatment <- round(N_patients / (N_treatments + 1))
  patients_no_treatment <- N_patients - patients_per_treatment * N_treatments
  
  

  #Put time effect into X
  if(time_effect) {
    N_fixed <- N_treatments + 1
    fixed_prior_sd <- array(
      c(array(c(prior$fixed_prior_sd_all), N_treatments), prior$time_effect_prior_sd),
      N_fixed)
  } else {
    N_fixed <- N_treatments
    fixed_prior_sd <- array(c(prior$fixed_prior_sd_all), N_treatments)
  }

  X <- array(0, c(N_patients, N_time, N_fixed))
  if(N_unknown_shift > 0) {
    n_shifts <- max_observation_shift + 1
  } else {
    n_shifts <- 0
  }
  X_unknown_shift <- array(0, c(N_unknown_shift, n_shifts, N_time, N_fixed))

  if(time_effect) {
    time_effect_shift <- -1
    for(p in 1:N_patients) {
      if(p %in% unknown_shift_patients) {
        for(time_shift_id in 1:max_observation_shift) {
          X_unknown_shift[unknown_shift_ids[p], time_shift_id, , N_fixed]  <- 1:N_time + time_effect_shift
        }
      } else {
        X[p,,N_fixed] <- 1:N_time + time_effect_shift
      }
    }
  }
    
  #Put treatment effect into X
  if(N_treatments > 0) {
    
    treatment_per_patient <- integer(N_patients)
    treatment_per_patient[1:patients_no_treatment] <- 0
    for(t in 1:N_treatments) {
      start <- (patients_no_treatment + (t - 1) * patients_per_treatment) + 1
      patients_for_treatment <- start:(start + patients_per_treatment - 1)
      treatment_per_patient[patients_for_treatment] <- t
      for(p in patients_for_treatment) {
        if(p %in% unknown_shift_patients) {
          for(time_shift in 0:max_observation_shift) {
            start_time_shifted <- treatment_start_times_observed[p] + time_shift
            X_unknown_shift[unknown_shift_ids[p], time_shift + 1, (start_time_shifted + 1):N_time, t]  <- 1
          }
        } else {
          X[p, (treatment_start_times[p] + 1):N_time, t] <- 1
        }
      }
    }
  } else {
    ## All zeroes
    treatment_per_patient <- integer(N_patients)
  }
  
  beta <- rnorm(N_fixed, 0, sd = fixed_prior_sd)
  
  states_true <- array(0L, c(N_patients, N_time))
  o_types_full <- array(0L, c(N_patients, N_time))
  viral_load_full <- array(NA_real_, c(N_patients, N_time))

  # print(X)
  # print(transition_thresholds)
  # print(beta)
  
  for(p in 1:N_patients) {  
    state <- rdunif(1, ill_states_shift + 1, ill_states_shift + N_ill_states)
    for(t in 1:N_time) {
      # Transition model - no transitions from healthy or severe states
      if(state != s_healthy && state != s_severe) { 
        linpred <- 0
        if(N_fixed > 0) {
          linpred <- linpred + sum(X[p,t,] * beta)
        }
        
        transition_probabilities <- ordered_logistic_probs(linpred + state_intercepts[state - ill_states_shift], transition_thresholds)
        
        if(abs(sum(transition_probabilities) - 1) > 1e-8) {
          print(transition_probabilities)
          stop("Probabilities don't sum to 1")
        }
        
        state <- sample(1:N_states, size = 1, prob = transition_probabilities) 
      }
      states_true[p, t] <- state
      
      # Observation model
      if(state == s_severe) {
        o_types_full[p, t] <- o_severe
      } else if(state == s_healthy) {
        if(runif(1) < specificity) {
          o_types_full[p, t] <- o_neg
        } else {
          o_types_full[p, t] <- o_pos
          viral_load_mean <- ill_mean_viral_load[rdunif(1, a = 1, b = N_ill_states)]
          repeat {
            viral_load_full[p, t] <- rnorm(1, viral_load_mean, sd = observation_sigma)
            if(viral_load_full[p, t] > 0) {
              break;
            }
          }
        }
      } else {
        if(runif(1) < 1 - sensitivity) {
          o_types_full[p, t] <- o_neg
        } else {
          o_types_full[p, t] <- o_pos
          viral_load_mean <- ill_mean_viral_load[state - ill_states_shift]
          repeat {
            viral_load_full[p, t] <- rnorm(1, viral_load_mean, sd = observation_sigma)
            if(viral_load_full[p, t] > 0) {
              break;
            }
          }
        }
      }
    }
  }
  
  o_types = array(NA_integer_, N_obs)
  viral_load = array(0, N_obs)
  viral_load_known = array(NA_integer_, N_obs)
  for(n in 1:N_obs) {
    if(patients[n] %in% unknown_shift_patients) {
      observation_shift <- unknown_shift[unknown_shift_ids[patients[n]]]
    } else {
      observation_shift <- 0
    }
    
    
    o_types[n] <- o_types_full[patients[n], times[n]]
    viral_load_known[n] <- patients[n] > N_patients_unknown_load & o_types[n] != o_neg & o_types[n] != o_severe
    if(viral_load_known[n]) {
      viral_load[n] <- viral_load_full[patients[n], times[n]]
    }
    
    times[n] <- times[n] - observation_shift
  }

  list(
    observed = c(prior, 
                 transform_predictors_to_unique(X, X_unknown_shift),
                 list(
      
      N_patients = N_patients,
      N_obs = N_obs,
      N_time = N_time,
      N_fixed = N_fixed,
      N_ill_states = N_ill_states,
      central_ill_state = central_ill_state,
      use_severe_state = use_severe_state,
      ill_mean_viral_load = ill_mean_viral_load,
      
      N_unknown_shift = N_unknown_shift,
      unknown_shift_patients = unknown_shift_patients,
      max_observation_shift = max_observation_shift,
      
      o_neg = o_neg,
      o_pos = o_pos,
      o_severe = o_severe,
      o_types = o_types,
      patients = patients,
      times = times,
      viral_load = viral_load,
      viral_load_known = viral_load_known,

      fixed_prior_sd = fixed_prior_sd,
      generate_predictions = 0,
      N_new_patient_predictions = 0,
      template_patients_for_new_prediction = integer(0),
      
      # Used only in plot functions
      N_treatments = N_treatments,
      treatment_per_patient = treatment_per_patient,
      treatment_start_times = treatment_start_times,
      treatment_start_times_observed = treatment_start_times_observed
     )),
    true = list(
      beta = beta,
      sensitivity = sensitivity,
      specificity = specificity,
      transition_thresholds = transition_thresholds,
      state_intercepts_high = state_intercepts_high,
      neg_state_intercepts_low = neg_state_intercepts_low,
      observation_sigma = observation_sigma,
      observation_shift_pred = unknown_shift,
      # used only in plot functions
      states_true = states_true
      #observation_model = observation_model
    )
  )
}

plot_sim_data_observed_hmm <- function(data) {
  
  severe_state_viral_load <- 35
  
  viral_load_jitter <- 3
  
  patient_data <- data.frame(patient = 1:data$observed$N_patients, 
                             treatment_start_time = data$observed$treatment_start_times_observed,
                             treatment = data$observed$treatment_per_patient) 
  
  if(data$observed$N_treatments > 0) {
    treatment_labels <- c("No", paste0("t_",1:data$observed$N_treatments))
  } else {
    treatment_labels <- "No"
  }
  
  observed_dataset <- data.frame(patient = data$observed$patients, 
                                 viral_load_known = data$observed$viral_load_known,
                                 viral_load_raw = data$observed$viral_load) %>%
    mutate(
      o_type = factor(data$observed$o_types, levels = c(data$observed$o_neg,data$observed$o_pos, data$observed$o_severe), labels = c("NEG","POS", "SEVERE")),
      viral_load = case_when(viral_load_known == 1 ~ viral_load_raw,
                              o_type == "NEG" ~ -viral_load_jitter - 1,
                              o_type == "POS" ~ severe_state_viral_load / 2,
                              o_type == "SEVERE" ~ severe_state_viral_load + viral_load_jitter + 1,
       ) + runif(n(), -viral_load_jitter, viral_load_jitter),
       type = interaction(o_type, factor(data$observed$viral_load_known, levels = c(0,1), labels = c("unknown","known")))
    )  %>%
    inner_join(patient_data) %>%
    mutate(
      time = data$observed$times,
      treated = treatment > 0 & time >= treatment_start_time,
      treated_group = factor(treated * treatment, levels = 0:data$observed$N_treatments, labels = treatment_labels))
  
  observed_dataset <- observed_dataset %>%
    mutate(time = time + runif(n(), -0.2, 0.2)) %>%
    arrange(time)
  
  
  observed_dataset %>% 
    ggplot(aes(x = time, y = viral_load, group = patient, shape = type)) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = severe_state_viral_load, color = "blue", linetype = "dashed") +
    geom_line(aes(color = treated_group)) + geom_point() +
    expand_limits(y = c(-1 - 2 * viral_load_jitter, severe_state_viral_load + 2 * viral_load_jitter + 1)) 
}


plot_sim_data_true_hmm <- function(data) {
  
  patient_data <- data.frame(patient = 1:data$observed$N_patients, 
                             treatment_start_time = data$observed$treatment_start_times,
                             treatment = data$observed$treatment_per_patient) 
  
  if(data$observed$N_treatments > 0) {
    treatment_labels <- c("No", paste0("t_",1:data$observed$N_treatments))
  } else {
    treatment_labels <- "No"
  }

  s_healthy <- 1
  s_severe <- data$observed$N_ill_states + 2
  severe_state_viral_load <- 35
  
  viral_load_jitter <- 3
  
  true_dataset <- data.frame(patient = 1:data$observed$N_patients) %>%
    crossing(data.frame(time = 1:data$observed$N_time)) %>%
    rowwise() %>%
    mutate(state_true = data$true$states_true[patient, time],
           viral_load = case_when(state_true == s_healthy ~ -viral_load_jitter - 1,
                                  state_true == s_severe ~ severe_state_viral_load + viral_load_jitter + 1,
                                      # The pmax is just to let the expresssion evaluate correctly when state_true = 1
                                  TRUE ~ data$observed$ill_mean_viral_load[pmax(state_true - 1, 1)]
           )
    ) %>%
    ungroup() %>% 
    inner_join(patient_data) %>%
    mutate(
      treated = treatment > 0 & time >= treatment_start_time,
      treated_group = factor(treated * treatment, levels = 0:data$observed$N_treatments, labels = treatment_labels))
  
  true_dataset <- true_dataset %>%
    mutate(time = time + runif(n(), -0.2, 0.2),
           viral_load = viral_load + runif(n(), -viral_load_jitter, viral_load_jitter)) %>%
    arrange(time)
  
  
  true_dataset %>% ggplot(aes(x = time, y = viral_load, group = patient)) + 
    geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = severe_state_viral_load, color = "blue", linetype = "dashed") +
    geom_line(aes(color = treated_group)) + geom_point() +
    expand_limits(y = c(-1 - 2 * viral_load_jitter, severe_state_viral_load + 2 * viral_load_jitter + 1)) 
}

