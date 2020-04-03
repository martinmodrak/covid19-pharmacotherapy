prior <- list(
  observation_model_prior = array(c(10,2,2,10), c(2,2)),
  
  p_healthy_prior_q025 = 0.01,
  p_healthy_prior_q975 = 0.5,
  
  viral_load_intercept_prior_sd = 5,
  viral_load_intercept_prior_mean = 0,
  viral_load_sigma_prior_sd = 5,
  
  initial_viral_load_mu = 20,
  initial_viral_load_sd = 5,
  kalman_total_noise_prior_q025 = 0.5,
  kalman_total_noise_prior_q975 = 10,
  
  fixed_prior_sd_all = 1,
  time_effect_prior_sd = 0.5
)

calculate_derived_prior_values <- function(prior) {
  normal_logit_match_error <- function(p, q025, q975) {
    mu <- p[1]
    sd <- p[2]
    (pnorm(gtools::logit(q025), mu, sd) - 0.025)^2 + (pnorm(gtools::logit(q975), mu, sd) - 0.975)^2
  }
  
  opt_res <- nlm(normal_logit_match_error, p = c(0,1), 
                 q025 = prior$p_healthy_prior_q025, q975 = prior$p_healthy_prior_q975)
  prior$state_intercept_prior_mean = opt_res$estimate[1]
  prior$state_intercept_prior_sd = opt_res$estimate[2]
  
  inv_gamma_match_error <- function(p, q025, q975) {
    alpha <- p[1]
    beta <- p[2]
    inv_gamma_cdf <- function(x, alpha, beta) { pgamma(beta / x, alpha, lower = FALSE) }
    (inv_gamma_cdf(q025, alpha, beta) - 0.025)^2 + (inv_gamma_cdf(q975, alpha, beta) - 0.975)^2
  }
  
  opt_res <- nlm(inv_gamma_match_error, p = c(1,1), q025 = prior$kalman_total_noise_prior_q025, 
                 q975 = prior$kalman_total_noise_prior_q025 )
  prior$kalman_total_noise_prior_alpha = opt_res$estimate[1]
  prior$kalman_total_noise_prior_beta = opt_res$estimate[2]
  
  prior
}

simulator <- function(N_patients, N_patients_unknown_load, N_treatments, N_obs_per_patient, N_time, prior, 
                      time_effect = FALSE) {
  prior <- calculate_derived_prior_values(prior)
  
  # Setup the HMM part
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

  state_intercept <- rnorm(1, prior$state_intercept_prior_mean, prior$state_intercept_prior_sd)
  
  # Setup the Kalman part
  
  total_noise <- 1 / rgamma(1, shape = prior$kalman_total_noise_prior_alpha, rate = prior$kalman_total_noise_prior_beta)

  process_noise_frac <- runif(1)
  process_noise_sd = sqrt( (total_noise ^ 2) * process_noise_frac);
  obs_noise_sd = sqrt( (total_noise ^ 2) * (1 - process_noise_frac));
  
  viral_load_intercept <- rnorm(1, prior$viral_load_intercept_prior_mean, prior$viral_load_intercept_prior_sd)
  
  
  # Setup shared elements
  
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

  if(time_effect) {
    time_effect_shift <- -1
    for(p in 1:N_patients) {
      X[p,,N_fixed] <- 1:N_time + time_effect_shift
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
        X[p, (treatment_start_times[p] + 1):N_time, t] <- 1
      }
    }
  } else {
    ## All zeroes
    treatment_per_patient <- integer(N_patients)
  }
  
  beta <- rnorm(N_fixed, 0, sd = fixed_prior_sd)
  
  o_types_true <- array(0L, c(N_patients, N_time))
  o_types_full <- array(0L, c(N_patients, N_time))
  viral_load_full <- array(NA_real_, c(N_patients, N_time))
  viral_load_true_full <- array(NA_real_, c(N_patients, N_time))
  
  for(p in 1:N_patients) {  
    state <- s_ill
    viral_load_true <- rnorm(1, prior$initial_viral_load_mu, prior$initial_viral_load_sd)
    for(t in 1:N_time) {
      linpred <- sum(X[p,t,] * beta)
      
      # Update HMM
      if(state == s_ill) {
        p_ill_to_healthy <- 1 / (1 + exp(-linpred - state_intercept))
        if(runif(1) < p_ill_to_healthy) {
          state <- s_healthy
        }
      }
      o_types_true[p, t] <- state
      o_types_full[p, t] <- sample(1:N_obs_types, size = 1, prob = observation_model[state,])
      
      # Update Kalman
      viral_load_linpred <- viral_load_intercept + linpred
      viral_load_true <- viral_load_true + viral_load_linpred + rnorm(1, 0, process_noise_sd)
      viral_load_true_full[p, t] <- viral_load_true
      viral_load_full[p, t] <- rnorm(1, viral_load_true, obs_noise_sd)
    }
  }
  
  o_types = array(NA_integer_, N_obs)
  viral_load = array(0, N_obs)
  viral_load_known = array(NA_integer_, N_obs)
  for(n in 1:N_obs) {
    o_types[n] <- o_types_full[patients[n], times[n]]
    viral_load_known[n] <- patients[n] > N_patients_unknown_load & o_types[n] != o_neg
    if(viral_load_known[n]) {
      viral_load[n] <- viral_load_full[patients[n], times[n]]
    }
  }

  list(
    observed = c(prior, list(
      N_patients = N_patients,
      N_obs = N_obs,
      N_time = N_time,
      N_treatments = N_treatments,
      patients = patients,
      times = times,
      o_neg = o_neg,
      o_pos = o_pos,
      o_types = o_types,
      viral_load = viral_load,
      viral_load_known = viral_load_known,
      N_fixed = N_fixed,
      X = X,
      fixed_prior_sd = fixed_prior_sd,
      treatment_per_patient = treatment_per_patient,
      treatment_start_times = treatment_start_times
    )),
    true = list(
      beta = beta,
      state_intercept = state_intercept,
      viral_load_intercept = viral_load_intercept,
      kalman_total_noise = total_noise,
      kalman_process_noise_frac = process_noise_frac,
      sensitivity = sensitivity,
      specificity = specificity,
      viral_load_true_full = viral_load_true_full,
      o_types_true = o_types_true
      #observation_model = observation_model
    )
  )
}

plot_sim_data_observed <- function(data) {
  
  patient_data <- data.frame(patient = 1:data$observed$N_patients, 
                             treatment_start_time = data$observed$treatment_start_times,
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
      o_type = factor(data$observed$o_types, levels = c(data$observed$o_neg,data$observed$o_pos), labels = c("NEG","POS")),
      viral_load = case_when(viral_load_known == 1 ~ viral_load_raw,
                              o_type == "NEG" ~ min(viral_load_raw) - 4 + runif(n(), 0, 2),
                              o_type == "POS" ~ mean(viral_load_raw) - 1  + runif(n(), 0, 2)#,
                              #data$observed$o_types == 3 ~ max(viral_load_raw) + 2
       ),
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
  
  
  observed_dataset %>% ggplot(aes(x = time, y = viral_load, group = patient, shape = type)) + geom_line(aes(color = treated_group)) + geom_point()
}


plot_sim_data_true <- function(data) {
  
  patient_data <- data.frame(patient = 1:data$observed$N_patients, 
                             treatment_start_time = data$observed$treatment_start_times,
                             treatment = data$observed$treatment_per_patient) 
  
  if(data$observed$N_treatments > 0) {
    treatment_labels <- c("No", paste0("t_",1:data$observed$N_treatments))
  } else {
    treatment_labels <- "No"
  }
  
  true_dataset <- data.frame(patient = 1:data$observed$N_patients) %>%
    crossing(data.frame(time = 1:data$observed$N_time)) %>%
    rowwise() %>%
    mutate(viral_load = data$true$viral_load_true_full[patient, time],
           o_type_raw = data$true$o_types_true[patient, time]) %>%
    ungroup() %>% 
    mutate(
      o_type = factor(o_type_raw, levels = c(data$observed$o_neg,data$observed$o_pos), labels = c("NEG","POS")),
      type = o_type
      #type = interaction(o_type, factor(data$observed$viral_load_known, levels = c(0,1), labels = c("unknown","known")))
    )  %>%
    inner_join(patient_data) %>%
    mutate(
      treated = treatment > 0 & time >= treatment_start_time,
      treated_group = factor(treated * treatment, levels = 0:data$observed$N_treatments, labels = treatment_labels))
  
  true_dataset <- true_dataset %>%
    mutate(time = time + runif(n(), -0.2, 0.2)) %>%
    arrange(time)
  
  
  true_dataset %>% ggplot(aes(x = time, y = viral_load, group = patient, shape = type)) + geom_line(aes(color = treated_group)) + geom_point()
}

