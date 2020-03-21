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
  real<lower=0> recovery;
  real dfirst;
  real firstval_raw;
  real<lower=0> sigma;
}


transformed parameters {
  real<lower=0> firstval = exp(firstval_raw * initial_disease_prior_logsd + initial_disease_prior_logmean);  
  real slope = exp(-first) * dfirst * (exp(first) + exp(recovery));
  real initial = firstval - slope * log1p_exp(first - recovery);
}

model {
  //vector[N] mu = initial + slope * log1p_exp(to_vector(times) - recovery);
  vector[N] mu;
  for(i in 1:N) {
    if(times[i] == first) {
      mu[i] = firstval;
    } else {
      mu[i] = initial + slope * log1p_exp(times[i] - recovery);
    }
  }
  y ~ normal(mu, sigma);
  
  sigma ~ normal(0,1);
  //initial_raw ~ normal(0, 1);
  //slope ~ normal(0, 1);
  firstval_raw ~ normal(0, 1);
  recovery ~ gamma(baseline_recovery_shape, baseline_recovery_shape / baseline_recovery_mean);
  dfirst ~ normal(0, 1);
  //change_at_first_if_slope_1_raw ~ normal(0, 1);
  //helper ~ normal(0, 20);
}

