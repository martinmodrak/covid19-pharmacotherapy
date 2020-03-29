simulator <- function(N_time) {
  N_obs_types <- 2
  o_neg <-  1
  o_pos <-  2
  
  N_states <- 2
  s_healthy <- 1
  s_ill <- 2
  
  observation_model <- matrix(0, nrow = N_states, ncol = N_obs_types)
  observation_model[, 1] <- runif(N_states)
  observation_model[, 2] <- 1 - observation_model[, 1]
  
  p_ill_to_healthy <- runif(1)
  
  state <- s_ill
  observations <- integer(N_time)
  for(t in 1:N_time) {
    observations[t] <- sample(1:N_obs_types, size = 1, prob = observation_model[state,])
    if(state == s_ill) {
      if(runif(1) < p_ill_to_healthy) {
        state <- s_healthy
      }
    }
  }
  
  list(
    observed = list(
      N_time = N_time,
      observations = observations
    ),
    true = list(
      p_ill_to_healthy = p_ill_to_healthy,
      observation_model = observation_model
    )
  )
}