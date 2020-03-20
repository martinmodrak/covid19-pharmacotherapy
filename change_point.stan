data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int observation_patients[N_obs];
  int observation_type[N_obs];
  int observation_time[N_obs];
  vector[N_obs] observations;
  
  //Priors
  real<lower=0> sigma_prior_sd;
  real initial_disease_prior_logmean;
  real<lower=0> initial_disease_prior_logsd;
  real<lower=0> baseline_recovery_mean_prior_mean;
  real<lower=0> baseline_recovery_mean_prior_shape;
  real<lower=0> baseline_recovery_shape_prior_sd;
  real t_high_prior_logmean;
  real<lower=0> t_high_prior_logsd;
  real<lower=0> baseline_recovery_mean;
  real<lower=0> baseline_recovery_shape;
  real<lower=0> baseline_slopes_mean_prior_sd;
  real<lower=0> baseline_slopes_sd_prior_sd;
}

parameters {
  vector[N_patients] initial_disease_raw;
  vector<lower=0>[N_patients] baseline_recovery;
  real baseline_slopes_mean;
  real<lower=0> baseline_slopes_sd;
  vector[N_patients] baseline_slopes_raw;
  real<lower=0> sigma_raw;
  //real t_high_raw;
  real<lower=0> t_high;
}

transformed parameters {
  vector<lower=0>[N_patients] initial_disease = exp(initial_disease_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);
  //real<lower=0> t_high = exp(t_high_raw * t_high_prior_logsd + t_high_prior_logmean) * mean(initial_disease);
  real<lower=0> sigma = sigma_raw * sigma_prior_sd;
  vector[N_patients] baseline_slopes = baseline_slopes_raw * baseline_slopes_sd + baseline_slopes_mean;
}

model {
  initial_disease_raw ~ normal(0, 1);
  sigma_raw ~ normal(0, 1);
  t_high ~ lognormal(t_high_prior_logmean, t_high_prior_logsd);
  
  // baseline_recovery_mean ~ gamma(baseline_recovery_mean_prior_shape, baseline_recovery_mean_prior_shape / baseline_recovery_mean_prior_mean);
  // baseline_recovery_shape ~ lognormal(0, baseline_recovery_shape_prior_sd);
  baseline_recovery ~ gamma(baseline_recovery_shape, baseline_recovery_shape / baseline_recovery_mean);
  
  //t_high_raw ~ normal(0, 1);
  baseline_slopes_sd ~ normal(0, baseline_slopes_sd_prior_sd);
  baseline_slopes_mean ~ normal(0, baseline_slopes_mean_prior_sd);
  baseline_slopes_raw ~ normal(0, 1);
  
  for(o in 1:N_obs) {
    int patient = observation_patients[o];
    real mu = initial_disease[patient] + baseline_slopes[patient] * log1p_exp(observation_time[patient] - baseline_recovery[patient]);
    if(observation_type[o] == 0) {
      target += normal_lpdf(observations[o] | mu, sigma);
    } else if (observation_type[o] == -1) {
      //Below detection limit, which is 0
      target += normal_lcdf(0 | mu, sigma);
    } else if (observation_type[o] == 1) {
      //Only know it is positive, e.g. above detection limit
      target += normal_lccdf(0 | mu, sigma);
    } else if (observation_type[o] == 2) {
      //Above t_high
      target += normal_lccdf(t_high | mu, sigma);
    } else {
      reject("Unrecognized observation type");
    }
  }
  
}

