hmm_model_defs_df <- tibble(model_type = "hmm", N_ill_states = c(3, 5)) %>%
  crossing(tibble(prior_type = c("default", "wide"))) %>%
  crossing(tibble(time_0 = c("symptom_onset_impute", "symptom_onset_estimate" , "first_measurement"))) %>%
  crossing(tibble(include_no_followup = c(TRUE, FALSE))) %>%
  crossing(tibble(no_followup_azithromycin = c(TRUE,FALSE))) %>%
  crossing(tibble(time_effect = c(TRUE, FALSE))) 
  

# Remove nonsensical combinations, build names
hmm_model_defs_df <- hmm_model_defs_df %>%
  filter(include_no_followup | !no_followup_azithromycin) %>%
  mutate(name = paste0(model_type, "_", N_ill_states, "_", prior_type, "_", time_0, 
                      if_else(include_no_followup, "", "_exclude_nofollowup"),
                      if_else(no_followup_azithromycin, "_AZ", ""),
                      if_else(time_effect, "_t",""))
                      )

if(length(unique(hmm_model_defs_df$name)) != nrow(hmm_model_defs_df)) {
  stop("Names are not unique")
}


hmm_model_defs <- list()
for(i in 1:nrow(hmm_model_defs_df)) {
  hmm_model_defs[[hmm_model_defs_df$name[i]]] <- hmm_model_defs_df[i,] %>% as.list()
}

all_model_defs <- hmm_model_defs