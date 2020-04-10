prepare_data_for_pp_checks_hmm <- function(fit_list) {
  data_for_model <- fit_list$data_for_model
  data_wide <- fit_list$data_wide
  fit <- fit_list$fit
  

  observed_dataset <- prepare_observed_dataset(data_for_model, data_wide) %>%
    mutate(observed_id = 1:n())
  
  
  # Keep only those fitted values that correspond to observed values
  fitted_dataset <- prepare_fitted_dataset_hmm(data_for_model, data_wide, fit)
  pp_check_dataset <- fitted_dataset %>%
    inner_join(observed_dataset, by = c("patient" = "patient", "time" = "time")) %>%
    ungroup()
  
  list(
    data_for_model = data_for_model,
    predicted = pp_check_dataset,
    observed = observed_dataset
  )
}

filter_pp_check_data_hmm <- function(.data, ...) {
  list(
    data_for_model = .data$data_for_model,
    predicted = .data$predicted %>% filter(...),
    observed = .data$observed %>% filter(...)
  )
}

mutate_pp_check_data_hmm <- function(.data, ...) {
  list(
    data_for_model = .data$data_for_model,
    predicted = .data$predicted %>% mutate(...),
    observed = .data$observed %>% mutate(...)
  )
}


pp_check_hmm <- function(pp_check_dataset, observed_column, predicted_column, n_draws = NULL, group = NULL, ...) {
  if(is.null(n_draws)) {
    predicted_draws <- pp_check_dataset$predicted
  } else {
    draws <- sample(unique(pp_check_dataset$predicted$.draw), size = n_draws)
    predicted_draws <- pp_check_dataset$predicted %>% filter(.draw %in% draws)
  }
  yrep_df <- predicted_draws %>% select(observed_id, .draw, {{predicted_column}}) %>%
    pivot_wider(id_cols = observed_id, names_from = .draw, names_prefix = "d", values_from = {{predicted_column}}) %>%
    arrange(observed_id)
  
  observed <- pp_check_dataset$observed %>% arrange(observed_id)
  
  if(!identical(as.integer(observed$observed_id), as.integer(yrep_df$observed_id))) {
    stop("Data do not match")
  }
  
  y <- observed %>% pull({{observed_column}})
  yrep <- yrep_df %>% select(-observed_id) %>% as.matrix() %>% t()
  if(is.null(group)) {
    bayesplot::pp_check(y, yrep, ...)
  } else {
    bayesplot::pp_check(y, yrep, group = observed[[group]], ...)
  }
}


pp_check_viral_load_hmm <- function(pp_check_dataset, ...) {
  pp_check_hmm(pp_check_dataset %>% filter_pp_check_data_hmm(viral_load_known), viral_load, viral_load_pred, ...)
  
}

pp_check_o_neg_hmm <- function(pp_check_dataset, ...) {
  prop_neg <- function(x) {
    mean(x == pp_check_dataset$data_for_model$o_neg)
  }
  pp_check_hmm(pp_check_dataset, o_types, o_type_pred, stat = prop_neg, ...)
  
}

pp_check_o_pos_hmm <- function(pp_check_dataset, ...) {
  prop_pos <- function(x) {
    mean(x == pp_check_dataset$data_for_model$o_pos)
  }
  pp_check_hmm(pp_check_dataset, o_types, o_type_pred, stat = prop_pos, ...)
  
}