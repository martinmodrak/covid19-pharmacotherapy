softclamp <- function(x, low, high) {
  if_else(x < (high + low) / 2,
          low + log(1 + exp(x - low)),
          high - log(1 + exp(-x + high)))
}

simulator <- function(N_patients, N_patients_binary, N_treatments, N_obs_per_patient, N_time, prior) {
  N_fails <- 0
  repeat {
    
    sigma <- abs(rnorm(1, 0, prior$sigma_prior_sd))
    # t_high_raw <- rlnorm(1, prior$t_high_prior_logmean, prior$t_high_prior_logsd)
    # t_high <- t_high_raw * mean(initial_disease)
    t_high <- rlnorm(1, prior$t_high_prior_logmean, prior$t_high_prior_logsd)
    
    # initial_disease <- rlnorm(N_patients, prior$initial_disease_prior_logmean, prior$initial_disease_prior_logsd)
    
    initial_disease <- rnorm(N_patients, prior$initial_disease_prior_mean, prior$initial_disease_prior_sd)
    
    N_obs <- N_patients * N_obs_per_patient
    observation_patients <- rep(1:N_patients, each = N_obs_per_patient)
    
    observation_times <- integer(N_obs)
    treatment_start_times <- integer(N_patients)
    for(p in 1:N_patients) {
      first_time <- rdunif(1, a = 1, b = N_time -  2 * N_obs_per_patient + 1)
      #first_time <- 1
      times <- sort(sample(first_time:(first_time + 2 * N_obs_per_patient - 1), N_obs_per_patient))
      observation_times[observation_patients == p] <- times
      
      treatment_start_times[p] = rdunif(1, a = first_time, b = max(times) - 1)
    }
    
    baseline_slopes_sd <- abs(rnorm(1, 0, prior$baseline_slopes_sd_prior_sd))
    baseline_slopes_mean <- rnorm(1, 0, prior$baseline_slopes_mean_prior_sd)
    baseline_slopes <- rnorm(N_patients, baseline_slopes_mean, baseline_slopes_sd)
    
    treatment_slopes <- rnorm(N_treatments, 0, prior$treatment_slopes_prior_sd)
    treatment_design_matrix <- matrix(0, nrow = N_obs, ncol = N_treatments)
    patients_per_treatment <- round(N_patients / (N_treatments + 1))
    patients_no_treatment <- N_patients - patients_per_treatment * N_treatments
    if(N_treatments > 0) {
      treatment_per_patient <- integer(N_patients)
      treatment_per_patient[1:patients_no_treatment] <- 0
      for(t in 1:N_treatments) {
        start <- (patients_no_treatment + (t - 1) * patients_per_treatment) + 1
        patients_for_treatment <- start:(start + patients_per_treatment - 1)
        treatment_per_patient[patients_for_treatment] <- t
        observation_indices <- observation_patients %in% patients_for_treatment
        treatment_design_matrix[observation_indices, t] <- 
          observation_times[observation_indices] - treatment_start_times[observation_patients[observation_indices]]
      }
    } else {
      ## All zeroes
      treatment_per_patient <- integer(N_patients)
    }
    #No effect of treatment before treatment start
    treatment_design_matrix[treatment_design_matrix < 0 ] <- 0
    
    observation_noise <- rnorm(length(observation_patients), mean = 0, sd = sigma)
    observation_true_value <- initial_disease[observation_patients] + baseline_slopes[observation_patients] * observation_times + 
      as.numeric(treatment_design_matrix %*% treatment_slopes) +
      observation_noise
    
    # clamp_shift = 5 * sigma + 2;
    # clamp_lower = -clamp_shift;
    # clamp_upper = t_high + clamp_shift;
    # observation_true_value <- softclamp(observation_true_value, clamp_lower, clamp_upper)
    
    
    observation_type <- case_when(observation_true_value < 0 ~ -1,
                                  observation_true_value > t_high ~ 2,
                                  observation_patients <= N_patients_binary ~ 1, #The POS case
                                  TRUE ~ 0)
    # observation_type = rep(0, length.out = N_obs)
    
    
    
    observations <- observation_true_value
    observations[observation_type != 0] <- 0
    
    res <- list(
      observed = c(prior, list(
        initial_disease = initial_disease,
        
        N_patients = N_patients,
        N_obs = N_obs,
        N_treatments = N_treatments,
        observation_patients = observation_patients,
        observation_time = observation_times,
        observation_type = observation_type,
        observations = observations,
        treatment_design_matrix = treatment_design_matrix,
        treatment_start_times = treatment_start_times,
        treatment_per_patient = treatment_per_patient
      )),
      true = list(
        sigma = sigma,
        t_high = t_high,
        baseline_slopes_sd = baseline_slopes_sd,
        baseline_slopes_mean = baseline_slopes_mean,
        baseline_slopes = baseline_slopes,
        treatment_slopes = treatment_slopes
      )
    )
    if(sum(observation_type == 0) >= 3) {
      break;
    }
    N_fails <- N_fails + 1
  }
  if(N_fails > 0){
    cat(N_fails, " fails\n")
  }
  res
}

plot_sim_data_true <- function(data) {
  
  patient_data <- data.frame(patient = 1:data$observed$N_patients, initial = data$observed$initial_disease, 
                             treatment_start_time = data$observed$treatment_start_times,
                             treatment = data$observed$treatment_per_patient,
                             slope = data$true$baseline_slopes) 
  
  treatment_slopes_shifted <- c(0, data$true$treatment_slopes)
  
  if(data$observed$N_treatments > 0) {
    treatment_labels <- c("No", paste0("t_",1:data$observed$N_treatments))
  } else {
    treatment_labels <- "No"
  }
  true_dataset <- patient_data %>% 
    crossing(time = seq(1, max(data$observed$observation_time) + 1, length.out = 50))  %>%
    mutate(
      treated = treatment > 0 & time > treatment_start_time,
      expression = initial + slope * time +
        if_else(treated,
                (time - treatment_start_time) * treatment_slopes_shifted[treatment + 1], 0),
      treated_group = factor(treated * treatment, levels = 0:data$observed$N_treatments, labels = treatment_labels)
    )
  
  
  true_dataset %>% ggplot(aes(x = time, y = expression, group = patient, color = treated_group)) + geom_line() + 
    geom_hline(color = "blue", linetype = "dashed", yintercept = 0) + geom_hline(color = "blue", linetype = "dashed", yintercept = data$true$t_high) 
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
  
  observed_dataset <- data.frame(patient = data$observed$observation_patients, 
                                 expression = case_when(data$observed$observation_type <= 0 ~ data$observed$observations,
                                                        data$observed$observation_type == 1 ~ mean(data$observed$observations),
                                                        data$observed$observation_type == 2 ~ max(data$observed$observations) + 2
                                 ),
                                 type = factor(data$observed$observation_type))  %>%
    inner_join(patient_data) %>%
    mutate(
      time = data$observed$observation_time,
      treated = treatment > 0 & time >= treatment_start_time,
      treated_group = factor(treated * treatment, levels = 0:data$observed$N_treatments, labels = treatment_labels))
  
  observed_dataset <- observed_dataset %>%
    mutate(time = time + runif(n(), -0.2, 0.2)) %>%
    arrange(time)
  
  
  observed_dataset %>% ggplot(aes(x = time, y = expression, group = patient, shape = type)) + geom_line(aes(color = treated_group)) + geom_point()
}
