model_color <- "#3993dd"
data_color <- "#91171f"
decoration_color <- "#ca5310"

prior_hmm <- list(
  fixed_prior_sd_all = 5,
  observation_sigma_prior_sd = 8,
  transition_thresholds_prior_sd = 5,
  state_intercept_prior_sd = 3,
  time_effect_prior_sd = 1#0.3
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
  maX_unknown_shift_shift <- dim(X_unknown_shift)[2]
    
  X_index <- array(-1, c(N_patients, N_time))
  X_unknown_shift_index <- array(-1, c(N_unknown_start, maX_unknown_shift_shift, N_time))
  
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
  
  for(u in 1:N_unknown_start) {
    for(shift in 1:maX_unknown_shift_shift) {
      for(t in 1:N_time) {
        str_index <- paste0(X_unknown_shift[u, shift,t,], collapse = "__x__")
        if(!is.null(X_dict[[str_index]])) {
          X_unknown_shift_index[u, shift, t] <- X_dict[[str_index]]$id
        } else {
          X_unknown_shift_index[u, shift, t] <- next_id
          X_dict[[str_index]] <- list(id = next_id, value = X_unknown_shift[u, shift,t,])
          next_id <- next_id + 1
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

prepare_data_for_plotting_hmm <- function(data_for_model, data_wide, fit) {
  pos_approx <- mean(data_for_model$ill_mean_viral_load[data_for_model$central_ill_state])
  
  sigma_summary <- summary(fit, pars = "observation_sigma")$summary 
  
  severe_approx <- data_for_model$ill_mean_viral_load[data_for_model$N_ill_states] + 2 * sigma_summary[,"97.5%"] + 5
  
  data_wide <- data_wide %>% mutate(patient_label = paste0(ID, " - ", case_when(Hydroxychloroquine == "No" ~ "Control",
                                                                                Azithromycin == "No" ~ "HCQ",
                                                                                TRUE ~ "HCQ + AZ")),
                                    treated = Hydroxychloroquine == "Yes" | Azithromycin == "Yes")

  o_pos <- data_for_model$o_pos
  o_neg <- data_for_model$o_neg
  o_severe <- data_for_model$o_severe
  
  s_severe <- data_for_model$N_ill_states + 2
  s_healthy <- 1
  # Subset the observed data
  observed_dataset <- data.frame(patient = data_for_model$patients,
                                 o_types = data_for_model$o_types,
                                 time = data_for_model$times,
                                 viral_load = data_for_model$viral_load,
                                 viral_load_known = data_for_model$viral_load_known) %>% 
    mutate(viral_load_imputed = case_when(
        viral_load_known ~ viral_load,
        o_types == o_pos ~ pos_approx,
        o_types == o_severe ~ severe_approx,
        o_types == o_neg ~ -2,
        TRUE ~ NA_real_
    ),
    type = 
      interaction(
        factor(o_types, levels = c(o_neg, o_pos, o_severe), labels = c("NEG","POS","ICU/DEATH")),
        factor(viral_load_known, levels = c(FALSE, TRUE), labels = c("unknown", "known")))
        )  %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    arrange(time) 

  observation_shifts <- spread_draws(fit, observation_shift_pred[unkown_shift_id]) %>%
    select(-.chain, -.iteration) %>% 
    mutate(patient = data_for_model$unknown_shift_patients[unkown_shift_id],
           observation_shift_pred = as.integer(observation_shift_pred))
  

  fitted_dataset <- spread_draws(fit, 
                                 state_pred[patient, time], o_type_pred[patient, time], 
                                 viral_load_pred[patient, time]) %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    ungroup() %>%
    left_join(observation_shifts, by = c("patient" = "patient", ".draw" = ".draw")) %>%
    ungroup() %>%
    mutate(viral_load_pred_imputed = case_when(
      state_pred == s_severe ~ severe_approx,
      state_pred == s_healthy ~ -2,
      TRUE ~ viral_load_pred
    )) %>% 
    mutate(time = if_else(!is.na(observation_shift_pred), time - observation_shift_pred, time))

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
  data_wide_new_patients <- data_wide %>% 
    filter(NumericID %in% data_for_model$template_patients_for_new_prediction) 
  
  
  new_patients_dataset <- spread_draws(fit, new_patient_state_pred[new_patient_id, time])  %>%
    mutate(patient_template = data_for_model$template_patients_for_new_prediction[new_patient_id]) %>%
    inner_join(data_wide_new_patients, by = c("patient_template" = "NumericID"))
  
  
  list(
    observed_dataset = observed_dataset,
    fitted_dataset = fitted_dataset,
    max_time = data_for_model$N_time,
    severe_approx = severe_approx,
    states_data = states_data,
    new_patients_dataset = new_patients_dataset
  )
  
}

plot_fitted_patients_hmm <- function(prepared_data, patient_ids, type = "fitted") {
  
  
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
    geom_line(color = data_color) + geom_point(color = data_color) + facet_wrap(~patient_label)
  
  
}


plot_fitted_new_patients_hmm <- function(prepared_data) {
  
    

    summarised_fitted <- prepared_data$new_patients_dataset %>%
      crossing(prepared_data$states_data) %>%
      group_by(patient_template, patient_label, time, state, low, high) %>%
      summarise(count = sum(new_patient_state_pred == state)) %>%
      group_by(patient_template, patient_label, time) %>%
      mutate(prob = count / sum(count)) %>%
      ungroup()
    
    severe_approx <- max(prepared_data$states_data$low)

    treatment_start_data <- prepared_data$new_patients_dataset %>% 
      ungroup() %>%
      select(patient_template, patient_label, treated, Days_From_Onset_Imputed) %>%
      distinct() %>% filter(treated) %>%
      mutate(treatment_start = Days_From_Onset_Imputed)
    
    #print(head(treatment_start_data))
    
    summarised_fitted %>% 
      ggplot(aes(xmin = time - 0.5, xmax = time + 0.5, ymin = low, ymax = high, fill = prob)) + 
      geom_rect() + scale_fill_gradient(low = "white", high = model_color, limits = c(0,1)) +
      geom_hline(yintercept = 0, color = decoration_color, linetype = "dashed") +
      geom_hline(yintercept = severe_approx, color = decoration_color, linetype = "dashed") +
      geom_vline(data = treatment_start_data, aes(xintercept = treatment_start), color = decoration_color, linetype = "dashed") +
      facet_wrap(~patient_label)
  
  
}
