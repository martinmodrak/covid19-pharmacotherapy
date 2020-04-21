model_color <- "#3993dd"
data_color <- "#91171f"
decoration_color <- "#ca5310"

prior_hmm <- list(
  fixed_prior_sd_all = 5,
  observation_sigma_prior_sd = 8,
  transition_thresholds_prior_sd = 5,
  state_intercept_prior_sd = 3,
  time_effect_prior_sd = 1
)

prior_wide_hmm <- list(
  fixed_prior_sd_all = 15,
  observation_sigma_prior_sd = 15,
  transition_thresholds_prior_sd = 15,
  state_intercept_prior_sd = 15,
  time_effect_prior_sd = 5
)


inv_logit <- function(a) {
  if (a < 0) {
    exp_a = exp(a);
    if (a < log(2.2e-16)) {
      exp_a;
    }
    exp_a / (1 + exp_a);
  }
  1/(1 + exp(-a));
}

ordered_logistic_probs <- function(lambda, c) {
  N = length(c) + 1
  res <- array(NA_real_, N)
  res[1] <- 1 - inv_logit(lambda - c[1])
  if(N > 2) {
    for(i in 2:(N-1)) {
      res[i] <- inv_logit(lambda - c[i - 1]) - inv_logit(lambda - c[i])
    }
  }
  res[N] <- inv_logit(lambda - c[N - 1])
  
  res
}

transform_predictors_to_unique <- function(X, X_unknown_shift) {
  if(typeof(X) != "double") {
    stop("has to be double")
  }
  if(length(dim(X)) != 3) {
    stop("Expecting 3 dims - patients, time, effect")
  }
  
  N_patients <- dim(X)[1]
  N_time <- dim(X)[2]
  N_fixed <- dim(X)[3]

  N_unknown_start <- dim(X_unknown_shift)[1]
  max_observation_shift_id <- dim(X_unknown_shift)[2]
    
  X_index <- array(-1, c(N_patients, N_time))
  X_unknown_shift_index <- array(-1, c(N_unknown_start, max_observation_shift_id, N_time))
  
  X_dict <- list()
  next_id <- 1
  for(p in 1:N_patients) {
    for(t in 1:N_time) {
      str_index <- paste0(X[p,t,], collapse = "__x__")
      if(!is.null(X_dict[[str_index]])) {
        X_index[p,t] <- X_dict[[str_index]]$id
      } else {
        X_index[p,t] <- next_id
        X_dict[[str_index]] <- list(id = next_id, value = X[p,t,])
        next_id <- next_id + 1
      }
    }
  }
  
  if(N_unknown_start > 0) {
    for(u in 1:N_unknown_start) {
      for(shift_id in 1:max_observation_shift_id) {
        for(t in 1:N_time) {
          str_index <- paste0(X_unknown_shift[u, shift_id,t,], collapse = "__x__")
          if(!is.null(X_dict[[str_index]])) {
            X_unknown_shift_index[u, shift_id, t] <- X_dict[[str_index]]$id
          } else {
            X_unknown_shift_index[u, shift_id, t] <- next_id
            X_dict[[str_index]] <- list(id = next_id, value = X_unknown_shift[u, shift_id,t,])
            next_id <- next_id + 1
          }
        }
      }
    }
  }
      
  
  N_predictor_sets <- next_id - 1
  X_new <- array(NA_real_, c(N_predictor_sets, N_fixed))
  for(n in names(X_dict)) {
    X_new[X_dict[[n]]$id,] <- X_dict[[n]]$value
  }
  
  if(any(is.na(X_new))) {
    stop("Should not be NA")
  }
  
  list(
    N_predictor_sets =N_predictor_sets,
    X = X_new,
    X_index = X_index,
    X_unknown_shift_index = X_unknown_shift_index
  )
}

gautret_to_long <- function(gautret) {
  gautret %>% 
    pivot_longer(matches("D[0-6]"), names_to = "Day", values_to = "PCR", names_prefix = "D") %>% 
    filter(!is.na(PCR)) %>%
    mutate(Day = as.integer(Day))
}

hmm_new_predictions_labels <- c("Control", "HCQ", "HCQ+AZ")
hmm_control_prediction_id <- 1
hmm_hcq_prediction_id <- 2
hmm_hcq_az_prediction_id <- 3

hmm_hcq_effect_id <- 1
hmm_az_effect_id <- 2
hmm_time_effect_id <- 3

create_data_for_model <- function(wide_for_model, prior_hmm, use_time_effect = FALSE, 
                                  time_0 = "symptom_onset_estimate", max_observation_shift = 5,
                                  N_ill_states = 3, N_time = NULL) {
  
  o_neg = 1
  o_pos = 2
  o_severe = 3
  
  
  long_for_model <- wide_for_model %>% gautret_to_long()
  
  unknown_start_wide <- wide_for_model %>% filter(is.na(Days_From_Onset))
  
  observation_type <- case_when(long_for_model$PCR == "NEG" ~ o_neg,
                                long_for_model$PCR == "POS" ~ o_pos,
                                long_for_model$PCR  %in% c("ICU","DEATH") ~ o_severe,
                                TRUE ~ o_pos
  )
  
  N_fixed <- 2 + use_time_effect
  
  central_ill_state <- ceiling(N_ill_states / 2)
  
  data_for_model <- c(prior_hmm, list(
    N_patients = nrow(wide_for_model),
    N_obs = nrow(long_for_model),
    N_fixed = N_fixed,
    N_ill_states = N_ill_states,
    central_ill_state = central_ill_state,

    o_neg = o_neg,
    o_pos = o_pos,
    o_severe = o_severe,
    use_severe_state = any(observation_type == o_severe),
    
    patients = long_for_model$NumericID,
    o_types = observation_type,
    viral_load_known = !(long_for_model$PCR %in% c("POS","NEG","ICU","DEATH")),
    fixed_prior_sd = rep(1, N_fixed),
    generate_predictions = 1,
    N_new_predictions = 3
  ))
  
  if(time_0 %in% c("symptom_onset_estimate", "symptom_onset_impute")) {
    data_for_model$times = long_for_model$Day + long_for_model$Days_From_Onset_Imputed + 1
  } else if(time_0 == "first_measurement") {
    data_for_model$times = long_for_model$Day + 1
  }
  
  data_for_model$viral_load = rep(0, nrow(long_for_model))
  data_for_model$viral_load[data_for_model$viral_load_known] <- PCR_Neg_Limit - 
    as.double(long_for_model$PCR[data_for_model$viral_load_known])
  
  # Separate the range from 0 to viral_load_high_bound into N_ill_states sectors and choose the middle
  state_viral_load_bounds <- quantile(data_for_model$viral_load[data_for_model$viral_load_known == 1], 
           probs = seq(0,1, length.out = N_ill_states + 1))
  ill_mean_viral_load <- 0.5 * (state_viral_load_bounds[1 : N_ill_states] + state_viral_load_bounds[2 : (N_ill_states + 1)])
  dim(ill_mean_viral_load) <- N_ill_states
  data_for_model$ill_mean_viral_load <- ill_mean_viral_load
  
  if(is.null(N_time)) {
    data_for_model$N_time <- max(data_for_model$times)
  } else {
    data_for_model$N_time <- N_time
  }
  
  X <- array(0, c(data_for_model$N_patients, data_for_model$N_time, N_fixed))
  for(p in 1:data_for_model$N_patients) {
    patient_data <- wide_for_model %>% filter(NumericID == p)
    if(time_0 == "first_measurement") {
      treatment_times <- 1 : data_for_model$N_time
    } else {
      treatment_times <- (patient_data$Days_From_Onset_Imputed + 1) : data_for_model$N_time
    }
    if(patient_data$Hydroxychloroquine == "Yes") {
      X[p, treatment_times, hmm_hcq_effect_id] <- 1
    }
    if(patient_data$Azithromycin == "Yes") {
      X[p, treatment_times, hmm_az_effect_id] <- 1
    }
    if(use_time_effect) {
      X[p, ,hmm_time_effect_id] <- 0:(data_for_model$N_time - 1)
    }
  }
  
  if(time_0 == "symptom_onset_estimate" && nrow(unknown_start_wide) > 0) {
    data_for_model$N_unknown_shift = nrow(unknown_start_wide)
    data_for_model$unknown_shift_patients = unknown_start_wide$NumericID
    data_for_model$max_observation_shift = max_observation_shift
    X_unknown_shift = array(0, c(nrow(unknown_start_wide), max_observation_shift + 1, data_for_model$N_time, data_for_model$N_fixed))
    for(up in 1:nrow(unknown_start_wide)) {
      p <- unknown_start_wide$NumericID[up]
      patient_data <- wide_for_model %>% filter(NumericID == p)
      for(shift in 0:max_observation_shift) {
        treatment_start <- patient_data$Days_From_Onset_Imputed + 1 + shift
        treatment_times <- treatment_start : data_for_model$N_time
        if(patient_data$Hydroxychloroquine == "Yes") {
          X_unknown_shift[up, shift, treatment_times, hmm_hcq_effect_id] <- 1
        }
        if(patient_data$Azithromycin == "Yes") {
          X_unknown_shift[up, shift, treatment_times, hmm_az_effect_id] <- 1
        }
        if(use_time_effect) {
          X_unknown_shift[up, shift, , hmm_time_effect_id] <- 0:(data_for_model$N_time - 1)
        }
        
      }
    }
  } else {
    data_for_model$N_unknown_shift = 0
    data_for_model$unknown_shift_patients = integer(0)
    data_for_model$max_observation_shift = 0
    X_unknown_shift = array(0, c(0, 0, data_for_model$N_time, data_for_model$N_fixed))
  }
  
  data_for_model$N_new_predictions <- 3
  X_new_prediction <- array(0, c(3, data_for_model$N_time, N_fixed))
  X_new_prediction[hmm_hcq_prediction_id, , hmm_hcq_effect_id] <- 1
  X_new_prediction[hmm_hcq_az_prediction_id, , hmm_hcq_effect_id] <- 1
  X_new_prediction[hmm_hcq_az_prediction_id, , hmm_az_effect_id] <- 1
  if(use_time_effect) {
    for(np in 1:data_for_model$N_new_predictions)
      X_new_prediction[np,,3] <- 0:(data_for_model$N_time - 1)
  }
  data_for_model$X_new_prediction <- X_new_prediction
  
  c(data_for_model, transform_predictors_to_unique(X, X_unknown_shift))
}

prepare_observed_dataset <- function(data_for_model, data_wide) {
  o_pos <- data_for_model$o_pos
  o_neg <- data_for_model$o_neg
  o_severe <- data_for_model$o_severe
  
  
  observed_dataset <- data.frame(patient = data_for_model$patients,
                                 o_types = data_for_model$o_types,
                                 time = data_for_model$times,
                                 viral_load = data_for_model$viral_load,
                                 viral_load_known = data_for_model$viral_load_known) %>% 
    mutate(type = 
      interaction(
        factor(o_types, levels = c(o_neg, o_pos, o_severe), labels = c("NEG","POS","ICU/DEATH")),
        factor(viral_load_known, levels = c(FALSE, TRUE), labels = c("unknown", "known")))
    )  %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    arrange(time) 
  
  observed_dataset
}

prepare_fitted_dataset_hmm <- function(data_for_model, data_wide, fit) {
  # Subset the observed data
  

  
  fitted_dataset <- spread_draws(fit, 
                                 state_pred[patient, time], o_type_pred[patient, time], 
                                 viral_load_pred[patient, time]) 
  
  
  if(data_for_model$N_unknown_shift > 0) {
    observation_shifts <- spread_draws(fit, observation_shift_pred[unkown_shift_id]) %>%
      select(-.chain, -.iteration) %>% 
      mutate(patient = data_for_model$unknown_shift_patients[unkown_shift_id],
             observation_shift_pred = as.integer(observation_shift_pred)) %>%
      ungroup()
    
    fitted_dataset <- fitted_dataset %>%
      left_join(observation_shifts, by = c("patient" = "patient", ".draw" = ".draw")) %>%
      ungroup() %>%
      mutate(time = if_else(!is.na(observation_shift_pred), time - observation_shift_pred, time))
  }
  
  fitted_dataset
}

prepare_data_for_analysis_hmm <- function(fit_list) {
  data_for_model <- fit_list$data_for_model
  data_wide <- fit_list$data_wide
  fit <- fit_list$fit
  
  pos_approx <- mean(data_for_model$ill_mean_viral_load[data_for_model$central_ill_state])
  
  sigma_summary <- summary(fit, pars = "observation_sigma")$summary 
  
  if(data_for_model$use_severe_state) {
    severe_approx <- data_for_model$ill_mean_viral_load[data_for_model$N_ill_states] + 2 * sigma_summary[,"97.5%"] + 5
  } else {
    severe_approx <- NULL
  }
  
  o_pos <- data_for_model$o_pos
  o_neg <- data_for_model$o_neg
  o_severe <- data_for_model$o_severe
  
  
  observed_dataset <- prepare_observed_dataset(data_for_model, data_wide) %>%
    mutate(viral_load_imputed = case_when(
      viral_load_known ~ viral_load,
      o_types == o_pos ~ pos_approx,
      o_types == o_severe ~ severe_approx,
      o_types == o_neg ~ -2,
      TRUE ~ NA_real_
    ))

  s_severe <- data_for_model$N_ill_states + 2
  s_healthy <- 1
  
  fitted_dataset <- prepare_fitted_dataset_hmm(data_for_model, data_wide, fit) %>%
    mutate(viral_load_pred_imputed = case_when(
      state_pred == s_severe ~ severe_approx,
      state_pred == s_healthy ~ -2,
      TRUE ~ viral_load_pred
    )) %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    ungroup()
    
    
  ## Info about states
  
  N_ill_states <- data_for_model$N_ill_states
  diff_ill_mean <- diff(data_for_model$ill_mean_viral_load)
  inner_boundaries <- data_for_model$ill_mean_viral_load[1:(N_ill_states - 1)] + 
    diff_ill_mean / 2
  
  state_boundaries <- c(-5, 0, inner_boundaries, 
                        inner_boundaries[length(inner_boundaries)] + diff_ill_mean[length(diff_ill_mean)] , 35)
                        
  states_data <- data.frame(state = 1:(N_ill_states + 2), 
                            low = state_boundaries[1:(N_ill_states + 2)], 
                            high = state_boundaries[2:(N_ill_states + 3)])
    
  ## Predictions for new patients
  info_new_pred_wide <- tibble(pred_id = 1:length(hmm_new_predictions_labels), label = hmm_new_predictions_labels)
  
  
  new_pred_dataset <- spread_draws(fit, state_pred_new[pred_id, time])  %>%
    inner_join(info_new_pred_wide, by = c("pred_id" = "pred_id"))

  new_logp_dataset <- spread_draws(fit, log_p_new[pred_id, time, state])  %>%
    inner_join(info_new_pred_wide, by = c("pred_id" = "pred_id"))
  
  
  list(
    data_for_model = data_for_model,
    observed_dataset = observed_dataset,
    fitted_dataset = fitted_dataset,
    max_time = data_for_model$N_time,
    severe_approx = severe_approx,
    states_data = states_data,
    new_pred_dataset = new_pred_dataset,
    new_logp_dataset = new_logp_dataset
  )
  
}

plot_fitted_patients_hmm <- function(prepared_data, patient_ids, type = "fitted", ncol = 2) {
  
  
  fitted_dataset <- prepared_data$fitted_dataset %>% filter(patient %in% patient_ids) 
  observed_dataset <- prepared_data$observed_dataset %>% filter(patient %in% patient_ids) 
  
  if(type == "fitted") {
    summarised_fitted <- fitted_dataset %>%
      crossing(prepared_data$states_data) %>%
      group_by(patient, patient_label, time, state, low, high) %>%
      summarise(count = sum(state_pred == state)) %>%
      group_by(patient, patient_label, time) %>%
      mutate(prob = count / sum(count)) %>%
      ungroup()
    
    model_geom1 <- geom_rect(data = summarised_fitted, 
                              aes(xmin = time - 0.5, xmax = time + 0.5, ymin = low, ymax = high,
                                  fill = prob), inherit.aes = FALSE)
                                                            
    model_geom2 <- scale_fill_gradient(low = "white", high = model_color, limits = c(0,1))
    severe_approx <- max(prepared_data$states_data$low)
    observed_dataset <- observed_dataset %>% mutate(
      viral_load_imputed = if_else(viral_load_imputed > severe_approx, 
                                   severe_approx + 5,
                                   viral_load_imputed))
      
  } else if(type == "predicted") {
    fitted_draws <- sample(1:max(prepared_data$fitted_dataset$.draw), 50)
    
    
    model_geom1 <- geom_line(data = fitted_dataset %>% filter(.draw %in% fitted_draws), 
                             aes(x = time, y = viral_load_pred_imputed, group = .draw), 
                             alpha = 0.3, color = model_color, inherit.aes = FALSE)
    model_geom2 <- NULL
    
    severe_approx <- prepared_data$severe_approx
  } else {
    stop("Invalid type")
  }
  
  observed_dataset %>% ggplot(aes(x = time, y = viral_load_imputed, shape = type, group = patient)) + 
    model_geom1 + model_geom2 +
    geom_hline(yintercept = 0, color = decoration_color, linetype = "dashed") +
    geom_hline(yintercept = severe_approx, color = decoration_color, linetype = "dashed") +
    geom_line(color = data_color) + geom_point(color = data_color) + facet_wrap(~patient_label, ncol = ncol) +
    scale_shape_manual(values = c("NEG.unknown" = 0, "POS.unknown" = 1,"POS.known" = 16,"ICU/DEATH.unknown" = 15))
  
  
}


plot_fitted_new_patients_hmm <- function(prepared_data, ncol = 2) {

  summarised_fitted <- prepared_data$new_logp_dataset %>%
    inner_join(prepared_data$states_data, by = c("state" = "state")) %>%
    group_by(pred_id, label, time, state, low, high) %>%
    summarise(prob = mean(exp(log_p_new))) %>%
    ungroup()
  
  if(prepared_data$data_for_model$use_severe_state) {
    severe_approx <- max(prepared_data$states_data$low)
    severe_geom <- geom_hline(yintercept = severe_approx, color = decoration_color, linetype = "dashed")
  } else {
    severe_geom <- NULL
  }
  
  
  summarised_fitted %>% 
    ggplot(aes(xmin = time - 0.5, xmax = time + 0.5, ymin = low, ymax = high, fill = prob)) + 
    geom_rect() + scale_fill_gradient(low = "white", high = model_color, limits = c(0,1)) +
    geom_hline(yintercept = 0, color = decoration_color, linetype = "dashed") +
    severe_geom +
    facet_wrap(~label, ncol = ncol)
  
  
}


# get_samples_severe_risk_hmm <- function(prepared_data, ) {
# }
