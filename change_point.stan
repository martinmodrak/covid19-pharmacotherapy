functions {
  real softclamp(real low, real high, real x) {
    if(x < (high + low) / 2) {
      return low + log1p_exp(x - low);
    } else {
      return high - log1p_exp(-x + high);
    }
  }
}

data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int<lower=0> N_treatments;
  int observation_patients[N_obs];
  int observation_type[N_obs];
  int observation_time[N_obs];
  vector[N_obs] observations;

  matrix[N_obs, N_treatments] treatment_design_matrix;

  //Priors
  real<lower=0> sigma_prior_sd;
  real initial_disease_prior_mean;
  real<lower=0> initial_disease_prior_sd;
  // real initial_disease_prior_logmean;
  // real<lower=0> initial_disease_prior_logsd;
  real t_high_prior_logmean;
  real<lower=0> t_high_prior_logsd;
  real<lower=0> baseline_slopes_mean_prior_sd;
  real<lower=0> baseline_slopes_sd_prior_sd;
  real<lower=0> treatment_slopes_prior_sd;
  
  vector[N_patients] initial_disease;
}

transformed data {
  int first_measurement[N_patients] = rep_array(max(observation_time), N_patients);
  for(n in 1:N_obs) {
    first_measurement[observation_patients[n]] = min(first_measurement[observation_patients[n]], observation_time[n]);
  }
}

parameters {
  //vector[N_patients] initial_disease_raw;
  real baseline_slopes_mean;
  vector[N_patients] baseline_slopes_raw;
  real<lower=0> baseline_slopes_sd;
  vector[N_treatments] treatment_slopes;  

  real<lower=0> sigma;
  real<lower=0> t_high;
}

transformed parameters {
  //vector<lower=0>[N_patients] initial_disease = exp(initial_disease_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);
  //vector[N_patients] initial_disease = initial_disease_raw * initial_disease_prior_sd + initial_disease_prior_mean;
  vector[N_patients] baseline_slopes = baseline_slopes_raw * baseline_slopes_sd + baseline_slopes_mean;
}

model {
  //Used for clamping value outside truncation bounds to increase numerical stability
  real clamp_shift = 5 * sigma + 2;
  real clamp_lower = -clamp_shift;
  real clamp_upper = t_high + clamp_shift;
  
  //Linear model of treatment
  vector[N_obs] treatment_value;
  if(N_treatments > 0) {
    treatment_value = treatment_design_matrix * treatment_slopes;
  } else {
    treatment_value = rep_vector(0, N_obs);
  }
  
  //initial_disease_raw ~ normal(0, 1);
  sigma ~ normal(0, sigma_prior_sd);
  t_high ~ lognormal(t_high_prior_logmean, t_high_prior_logsd);
  
  //t_high_raw ~ normal(0, 1);
  baseline_slopes_sd ~ normal(0, baseline_slopes_sd_prior_sd);
  baseline_slopes_mean ~ normal(0, baseline_slopes_mean_prior_sd);
  baseline_slopes_raw ~ normal(0, 1);
  
  treatment_slopes ~ normal(0, treatment_slopes_prior_sd);
  
  for(n in 1:N_obs) {
    int patient = observation_patients[n];
    real mu = initial_disease[patient] + baseline_slopes[patient] * observation_time[n] + treatment_value[n];
    if(observation_type[n] == 0) {
      target += normal_lpdf(observations[n] | mu, sigma);
    } else {
      //real mu_clamp = softclamp(clamp_lower, clamp_upper, mu);
      real mu_clamp = mu;
      if (observation_type[n] == -1) {
        //Below detection limit, which is 0
        target += normal_lcdf(0 | mu_clamp, sigma);
      } else if (observation_type[n] == 1) {
        //Only know it is positive, e.g. above detection limit, also has to be below t_high
        target += log_diff_exp(normal_lcdf(t_high | mu_clamp, sigma), normal_lcdf(0 | mu_clamp, sigma));
      } else if (observation_type[n] == 2) {
        //Above t_high
        target += normal_lccdf(t_high | mu_clamp, sigma);
      } else {
        reject("Unrecognized observation type");
      }
    }
  }
}
