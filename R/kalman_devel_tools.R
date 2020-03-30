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
  
  #observed = rbinom(N_time, size = 1, p = 0.4)
  observed = rep(1,N_time)
  times = which(observed == 1);
  viral_load <- viral_load[observed == 1];
  N_obs = sum(observed);
  
  list(true = list(
    linpred = linpred,
    total_noise = total_noise,
    process_noise_frac = process_noise_frac
  ), observed = list(
    N_obs = N_obs,
    viral_load = viral_load,
    times = times,
    initial_viral_load_mu = initial_viral_load_mu,
    initial_viral_load_sd = initial_viral_load_sd
  ))
}