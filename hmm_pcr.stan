data {
  int<lower=0> N_time;
  //States: 1 = NEG, 2 = POS, 3 = ICU/DEATH
  int<lower=1,upper=2> observations[N_time];
}

transformed data {
  int N_obs_types = 2;
  int o_neg = 1;
  int o_pos = 2;
  int N_states = 2;
  int s_healthy = 1;
  int s_ill = 2;
  
}

parameters {
  real<lower=0,upper=1> p_ill_to_healthy;
  simplex[N_obs_types] observation_model[N_states];
}

model {
  //Assumming ill at time = 0
  matrix[N_time, N_obs_types] gamma;
  real acc_healthy[2];
  
  gamma[1, ] = transpose(observation_model[s_ill, ]);
  for(t in 2:N_time) {
    acc_healthy[s_healthy] = gamma[t - 1, s_healthy] + log(observation_model[s_healthy, observations[t]]);
    acc_healthy[s_ill] = gamma[t - 1, s_ill] + log(p_ill_to_healthy) + log(observation_model[s_healthy, observations[t]]);
    gamma[t, s_healthy] = log_sum_exp(acc_healthy);
    
    gamma[t, s_ill] =  gamma[t - 1, s_ill] + 
      log1m(p_ill_to_healthy) + 
      log(observation_model[s_ill, observations[t]]);
  }
  
  target += log_sum_exp(gamma[N_time,]);
}

