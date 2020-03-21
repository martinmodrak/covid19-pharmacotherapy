data {
  int<lower=0> N;
  vector[N] y;
  int times[N];
  real initial_disease_prior_logmean;
  real<lower=0> initial_disease_prior_logsd;
  real<lower=0> baseline_recovery_mean;
  real<lower=0> baseline_recovery_shape;
}

transformed data {
  int last = max(times);
  int first = min(times);
}

parameters {
  //real initial_raw;
  real change_at_first_if_slope_1_raw;
  //real helper;
  real slope;
  real lastval_raw;
  real<lower=0> sigma;
}


transformed parameters {
  //real<lower=0> initial = exp(initial_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);  
  real<lower=0> lastval = exp(lastval_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);  
  real<lower=0> change_at_first_if_slope_1 = exp(change_at_first_if_slope_1_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);  
  real helper = log_diff_exp(change_at_first_if_slope_1, 0);
  real recovery = first - helper;
  real initial = lastval - slope * log1p_exp(last - first + helper);
}

model {
  //vector[N] mu = initial + slope * log1p_exp(to_vector(times) - recovery);
  vector[N] mu;
  for(i in 1:N) {
    if(times[i] == last) {
      mu[i] = lastval;
    } else {
      mu[i] = initial + slope * (log1p_exp(times[i] - first + helper) - log1p_exp(last - first + helper));
    }
  }
  y ~ normal(mu, sigma);
  
  sigma ~ normal(0,1);
  //initial_raw ~ normal(0, 1);
  slope ~ normal(0, 1);
  lastval_raw ~ normal(0, 1);
  //recovery ~ gamma(baseline_recovery_shape, baseline_recovery_shape / baseline_recovery_mean);
  change_at_first_if_slope_1_raw ~ normal(0, 1);
  //helper ~ normal(0, 20);
}

