model_color <- "blue"

clamp_expression <- function(x, t_high) {
  case_when(x < -2 ~ runif(length(x), min = -3, max = -2),
            x > t_high ~ t_high,
            TRUE ~ x)
}

prepare_data_for_plotting <- function(data_for_model, data_wide, fit) {
  t_high_summary <- summary(fit, pars = "t_high")$summary 
  t_high_approx <- t_high_summary[,"50%"]
  pos_approx <- mean(data_for_model$observations)
  
  data_wide <- data_wide %>% mutate(patient_label = paste0(ID, " - ", case_when(Hydroxychloroquine == "No" ~ "Control",
                                                                                Azithromycin == "No" ~ "HCQ",
                                                                                TRUE ~ "HCQ + AZ")))
  
  # Subset the observed data
  observed_dataset <- data.frame(patient = data_for_model$observation_patients,
                                 observation_type = data_for_model$observation_type,
                                 time = data_for_model$observation_time,
                                 observation = data_for_model$observations) %>% 
    mutate(expression = case_when(observation_type <= 0 ~ observation,
                                  observation_type == 1 ~ pos_approx,
                                  observation_type == 2 ~ t_high_approx + 2
    ),
    type = factor(observation_type, levels = c(-1, 0, 1, 2), labels = c("NEG","Normal","POS","ICU/DEATH")))  %>%
    #filter(patient %in% patient_ids) %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    arrange(time) 
    
  max_time <- max(data_for_model$observation_time) + 2
  
  treatment_samples <- spread_draws(fit, treatment_slopes[treatment]) %>%
    pivot_wider(names_from = treatment, names_prefix = "treatment_", values_from = "treatment_slopes")
  
  fitted_dataset <- spread_draws(fit, initial_disease[patient], baseline_slopes[patient], t_high, sigma) %>% 
    #filter(patient %in% patient_ids) %>%
    inner_join(treatment_samples, by = c(".chain", ".iteration", ".draw")) %>%
    inner_join(data_wide, by = c("patient" = "NumericID")) %>%
    crossing(time = 0:max_time) %>%
    mutate(
      treatment_slope = case_when(time <= Days_From_Onset_Imputed ~ 0,
                                  Azithromycin == "Yes" ~ treatment_1 + treatment_2,
                                  Hydroxychloroquine == "Yes" ~ treatment_1,
                                  TRUE ~ 0),
      mu = initial_disease + time * baseline_slopes + 
        treatment_slope * (time - Days_From_Onset_Imputed),
      mu_clamped = clamp_expression(mu, t_high),
      y = rnorm(n(), mu, sigma),
      y_clamped = clamp_expression(y, t_high),
      treated = treatment_slope >= 0)
  
  list(
    observed_dataset = observed_dataset,
    fitted_dataset = fitted_dataset,
    max_time = max_time,
    t_high_summary = t_high_summary
    )
}

plot_fitted_patients <- function(prepared_data, patient_ids, type = "fitted") {

  fitted_draws <- sample(1:max(prepared_data$fitted_dataset$.draw), 50)
  
  fitted_dataset <- prepared_data$fitted_dataset %>% filter(patient %in% patient_ids, .draw %in% fitted_draws) 
  observed_dataset <- prepared_data$observed_dataset %>% filter(patient %in% patient_ids) 
  
  t_high_summary <- prepared_data$t_high_summary
  
  if(type == "fitted") {
    model_geom1 <- geom_line(data = fitted_dataset, aes(x = time, y = mu_clamped, color = treated, group = .draw), alpha = 0.3, inherit.aes = FALSE)
    model_geom2 <- NULL
  } else if(type == "predicted") {
    predicted_dataset <- fitted_dataset %>%
      group_by(patient, time) %>%
      summarise(low = quantile(y_clamped, probs = 0.025), low50 = quantile(y_clamped, probs = 0.25),
                high50 = quantile(y_clamped, probs = 0.75), high = quantile(y_clamped, probs = 0.975))
    
    
    model_geom1 <- geom_ribbon(data = predicted_dataset, aes(x = time, ymin = low,ymax = high, group = patient), fill = model_color, alpha = 0.3, inherit.aes = FALSE)
    model_geom2 <- geom_ribbon(data = predicted_dataset, aes(x = time, ymin = low50,ymax = high50, group = patient), fill = model_color, alpha = 0.3, inherit.aes = FALSE)
  } else {
    stop("Invalid type")
  }
  
  observed_dataset %>% ggplot(aes(x = time, y = expression, shape = type, group = patient)) + 
    model_geom1 + model_geom2 +
    geom_hline(yintercept = 0, color = "green", linetype = "dashed") +
    annotate("rect", xmin = 0, xmax = prepared_data$max_time, ymin = t_high_summary[,"2.5%"], ymax = t_high_summary[,"97.5%"], fill = "green", alpha = 0.3) +
    geom_line() + geom_point() + facet_wrap(~patient_label)
  
  
}