rstan_fit_with_cache <- function(cache_file, model, data, adapt_delta = 0.8) {
  if(inherits(model, "stanmodel")) {
    model_code <- model@model_code
  } else {
    model_code <- paste0(model$code(), collapse = "\n")
  }
  
  run_model <- TRUE
  if(file.exists(cache_file)) {
    cache_data <- readRDS(cache_file)
    if(is.list(cache_data) && !is.null(cache_data$model_code) 
       && !is.null(cache_data$data) && !is.null(cache_data$fit)) {
      normalized_code <- gsub("[[:space:]]", "", model_code)
      normalized_code_cache <- gsub("[[:space:]]", "", cache_data$model_code)
      if(identical(normalized_code, normalized_code_cache) && identical(data, cache_data$data)) {
        fit <- cache_data$fit
        run_model <- FALSE
      } else {
        message("Cache file exists, but is out of date")
      }
    } else {
      message("Cache file exists, but doesn't have usable contents")
    }
  }
  
  if(run_model) {
    if(inherits(model, "stanmodel")) {
      fit <- sampling(model, data = data, control = list(adapt_delta = adapt_delta))
      check_hmc_diagnostics(fit)
    } else {
      cmdstan_fit <- model$sample(data = data_for_model_time, adapt_delta = adapt_delta)
      fit <- rstan::read_stan_csv(cmdstan_fit$output_files())
    }
    
    saveRDS(list(model_code = model_code, data = data, fit = fit), cache_file)
  }
  
  fit
}

prepare_data_for_hmm_model <- function(data_wide, model_def) {
  if(model_def$model_type != "hmm") {
    stop("Wrong model def")
  }
  if(model_def$prior_type == "default") {
    prior <- prior_hmm
  } else if(model_def$prior_type == "wide") {
    prior <- prior_wide_hmm
  } else {
    stop("Invalid prior def")
  }

  if(!model_def$include_no_followup) {
    data_wide <- data_wide %>% filter(Source != "NoFollowup")
  }
  
  if(model_def$no_followup_azithromycin) {
    data_wide <- data_wide %>% mutate(Azithromycin = if_else(Source == "NoFollowup", "Yes", Azithromycin))
  }
  
  data_wide <- data_wide %>% mutate(NumericID = 1:n())

  
  list(data_wide = data_wide,
       data_for_model = 
        create_data_for_model(
          wide_for_model = data_wide, 
          prior_hmm = prior, 
          time_0 = model_def$time_0,
          use_time_effect = model_def$time_effect,
          max_observation_shift = 5,
          N_ill_states = model_def$N_ill_states
        )
  )
}

fit_hmm_model <- function(hmm_stanmodel, data_wide, model_def) {
  if(model_def$model_type != "hmm") {
    stop("Wrong model def")
  }
  
  data_list <- prepare_data_for_hmm_model(data_wide, model_def)
  
  cache_dir <- here("stored_fits")
  if(!dir.exists(cache_dir)) {
    dir.create(cache_dir)
  }
  
  fit <- rstan_fit_with_cache(paste0(cache_dir, "/", model_def$name, ".rds"), hmm_stanmodel, data_list$data_for_model, adapt_delta = 0.95)
  
  c(list(fit = fit), data_list)
}

fit_model <- function(hmm_stanmodel, data_wide, model_def) {
  if(model_def$model_type == "hmm") {
    fit_hmm_model(hmm_stanmodel, data_wide, model_def)
  } else {
    stop("Unknown model type")
  }
}

fit_all_models <- function(hmm_stanmodel, data_wide) {
  all_model_defs %>% purrr::map(fit_model, hmm_stanmodel = hmm_stanmodel, data_wide = data_wide)
}
