data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int observation_patients[N_obs];
  int observation_type[N_obs];
  int observation_time[N_obs];
  vector[N_obs] observations;
  int<lower=1> max_time_for_predictions;
  
  //Priors
  real<lower=0> sigma_prior_sd;
  real initial_disease_prior_logmean;
  real<lower=0> initial_disease_prior_logsd;
  real t_high_prior_logmean;
  real<lower=0> t_high_prior_logsd;
  real<lower=0> baseline_slopes_mean_prior_sd;
  real<lower=0> baseline_slopes_sd_prior_sd;
}

parameters {
  vector[N_patients] initial_disease_raw;
  real baseline_slopes_mean;
  vector[N_patients] baseline_slopes_raw;
  real<lower=0> baseline_slopes_sd;
  real<lower=0> sigma;
  //real t_high_raw;
  real<lower=0> t_high;
}

transformed parameters {
  vector<lower=0>[N_patients] initial_disease = exp(initial_disease_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);
  //real<lower=0> t_high = exp(t_high_raw * t_high_prior_logsd + t_high_prior_logmean) * mean(initial_disease);
  vector[N_patients] baseline_slopes = baseline_slopes_raw * baseline_slopes_sd + baseline_slopes_mean;
  
}

model {
  //Used for clamping value outside truncation bounds to increase numerical stability
  real clamp_shift = 5 * sigma + 2;
  real clamp_lower = -clamp_shift;
  real clamp_upper = t_high + clamp_shift;
  
  
  initial_disease_raw ~ normal(0, 1);
  sigma ~ normal(0, sigma_prior_sd);
  t_high ~ lognormal(t_high_prior_logmean, t_high_prior_logsd);
  
  //t_high_raw ~ normal(0, 1);
  baseline_slopes_sd ~ normal(0, baseline_slopes_sd_prior_sd);
  baseline_slopes_mean ~ normal(0, baseline_slopes_mean_prior_sd);
  baseline_slopes_raw ~ normal(0, 1);
  
  for(n in 1:N_obs) {
    int patient = observation_patients[n];
    real mu = initial_disease[patient] + baseline_slopes[patient] * observation_time[n];
    if(observation_type[n] == 0) {
      target += normal_lpdf(observations[n] | mu, sigma);
    } else {
      real mu_clamp_lower = clamp_lower + log1p_exp(mu - clamp_lower);
      real mu_clamp_upper = clamp_upper - log1p_exp(-mu + clamp_upper);
      
      if (observation_type[n] == -1) {
        //Below detection limit, which is 0
        target += normal_lcdf(0 | mu_clamp_lower, sigma);
      } else if (observation_type[n] == 1) {
        //Only know it is positive, e.g. above detection limit, also has to be below t_high
        target += log_diff_exp(normal_lcdf(t_high | mu_clamp_upper, sigma), normal_lcdf(0 | mu_clamp_lower, sigma));
      } else if (observation_type[n] == 2) {
        //Above t_high
        target += normal_lccdf(t_high | mu_clamp_upper, sigma);
      } else {
        reject("Unrecognized observation type");
      }
    }
  }
  
}

generated quantities {
  matrix[N_patients, max_time_for_predictions] mu_pred;
  matrix[N_patients, max_time_for_predictions] y_pred;
  
  for(patient in 1:N_patients) {
    for(t in 1:max_time_for_predictions) {
      real mu = initial_disease[patient] + baseline_slopes[patient] * t;
      mu_pred[patient, t] = mu;
      y_pred[patient, t] = normal_rng(mu, sigma);
    }
  }
}
