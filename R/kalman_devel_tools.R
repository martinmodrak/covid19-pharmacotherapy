simulator_kalman <- function(N_time) {
  initial_viral_load_mu <- 10
  initial_viral_load_sd <- 5
  
  linpred <- rnorm(1)
  
  process_noise_frac <- runif(1)
  total_noise <- abs(rnorm(1))
  process_noise = sqrt( (total_noise ^ 2) * process_noise_frac);
  sigma = sqrt( (total_noise ^ 2) * (1 - process_noise_frac));
  
  # sigma <- abs(rnorm(1))
  # process_noise <- abs(rnorm(1))

  state <- rnorm(1, initial_viral_load_mu, initial_viral_load_sd)
  viral_load <- numeric(N_time)
  for(t in 1:N_time) {
    state <- state + linpred + rnorm(1, 0, process_noise)
    viral_load[t] <- rnorm(1, state, sigma)
  }
  
  list(true = list(
    linpred = linpred,
    total_noise = total_noise,
    process_noise_frac = process_noise_frac
  ), observed = list(
    N_time = N_time,
    viral_load = viral_load,
    initial_viral_load_mu = initial_viral_load_mu,
    initial_viral_load_sd = initial_viral_load_sd
  ))
}