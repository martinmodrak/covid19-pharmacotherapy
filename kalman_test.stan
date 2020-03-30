data {
  int<lower=0> N_obs;
  int times[N_obs];
  vector[N_obs] viral_load;
  real initial_viral_load_mu;
  real initial_viral_load_sd;
}

transformed data {
  int step_size[N_obs];
  step_size[1] = times[1];
  for(n in 2:N_obs){
    if(times[n] <= times[n - 1]) {
      reject("Time must increase");
    }
    step_size[n] = times[n] - times[n - 1];
  }
}

parameters {
  real<lower=0> total_noise;
  real<lower=0,upper=1> process_noise_frac;
  real linpred;
}

transformed parameters {
  real<lower=0> process_variance = (total_noise ^ 2) * process_noise_frac;
  real<lower=0> obs_variance = (total_noise ^ 2) * (1 - process_noise_frac);

}

model {

  vector[N_obs] viral_load_centered;
  for(t in 1:N_obs) {
    viral_load_centered[t] = viral_load[t] - (times[t] - 1) * linpred;
  }
  
  {
    //Code based on the implementation in
    // https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp
    // Using the same notation, simplifed because both design matrix (F) and transition matrix (G) are
    // identity
    
    real m = initial_viral_load_mu;
    real C = initial_viral_load_sd ^ 2;
    for(n in 1:N_obs) {
      real Q;
      real Q_inv;
      real e;
      real A;
      C += process_variance * step_size[n];
      //f = m;
      Q = C + obs_variance;
      Q_inv = inv(Q);
      A = C * Q_inv;
      //C -= Q * A * A;
      C *= 1 - A;
      e = viral_load_centered[n] - m;
      m += A*e;
      target += -0.5 * (log(Q) + e * e * Q_inv);
    }
  }
  
  // target += gaussian_dlm_obs_lpdf(to_matrix(viral_load_centered, 1, N_time)| to_matrix({1}, 1, 1), to_matrix({1}, 1, 1),
  //   to_matrix({obs_variance}, 1, 1), to_matrix(process_variance}, 1,1), to_vector({initial_viral_load_mu}),
  //   to_matrix({initial_viral_load_sd ^ 2}, 1, 1));
  
  linpred ~ normal(0, 1);
  total_noise ~ normal(0,1);
  // sigma ~ normal(0, 1);
  // process_noise ~ normal(0, 1);
}

