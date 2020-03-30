data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int<lower=0> N_fixed;
  
  //Observations: 1 = NEG, 2 = POS, 3 = ICU/DEATH
  int<lower=1, upper=2> o_types[N_obs];
  int<lower=1, upper=N_patients> patients[N_obs];
  int<lower=1> times[N_obs];
  vector[N_obs] viral_load;
  matrix[N_obs, N_fixed] X;
  
  
  real state_intercept_prior_sd;
  real viral_load_intercept_prior_sd;
  real viral_load_intercept_prior_mean;
  real viral_load_sigma_prior_sd;
  
  vector<lower=0>[N_fixed] fixed_prior_sd;
  real time_effect_prior_sd;
  real time_effect_shift;
}

transformed data {
  int N_obs_types = 2;
  int o_neg = 1;
  int o_pos = 2;
  int N_states = 2;
  int s_healthy = 1;
  int s_ill = 2;
  int N_time = max(times);
  //0 means no observation
  int<lower=0,upper=N_obs> obs_ids[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0> patient_max_time[N_patients];
  
  for(n in 1:N_obs) {
    patient_max_time[patients[n]] = max(times[n], patient_max_time[patients[n]]);
    obs_ids[patients[n], times[n]] = n;
  }
}

parameters {
  vector[N_fixed] beta;
  real state_intercept;
  real viral_load_intercept;
  real<lower=0> viral_load_sigma;
  real time_effect;
  real<lower=0.5, upper=1> sensitivity; 
  real<lower=0.5, upper=1> specificity; 
}

transformed parameters {
  simplex[N_obs_types] observation_model[N_states];
  observation_model[s_healthy, o_neg] = specificity; 
  observation_model[s_healthy, o_pos] = 1 - specificity;
  observation_model[s_ill, o_neg] = 1 - sensitivity; 
  observation_model[s_ill, o_pos] = sensitivity;
}

model {
  //Assumming ill at time = 0
  matrix[N_time, N_obs_types] log_forward_p;
  //matrix[N_time, N_obs_types] log_backward_p;
  real acc_healthy[2];
  real acc_ill[2];
  
  //log_backward_p[N_time, ] = rep_row_vector(0, N_obs_types);
  for(p in 1:N_patients) {
    if(obs_ids[p, 1] == 0) {
      log_forward_p[1, s_ill] = 0;
    } else {
      log_forward_p[1, s_ill] = log(observation_model[s_ill, o_types[obs_ids[p, 1]]]);
    }
    log_forward_p[1, s_healthy] = log(0);
    for(t in 2:patient_max_time[p]) {
      int id = obs_ids[p, t];
      real base_linpred = X[id,] * beta + time_effect * (t + time_effect_shift);
      real p_ill_to_healthy = inv_logit(state_intercept + base_linpred);
      acc_healthy[s_healthy] = log_forward_p[t - 1, s_healthy];
      acc_healthy[s_ill] = log_forward_p[t - 1, s_ill] + log(p_ill_to_healthy);
      log_forward_p[t, s_healthy] = log_sum_exp(acc_healthy);
      if(id != 0) {
        log_forward_p[t, s_healthy] += log(observation_model[s_healthy, o_types[id]]);
      }
      
      log_forward_p[t, s_ill] =  log_forward_p[t - 1, s_ill] + 
        log1m(p_ill_to_healthy);
      if(id != 0) {
        log_forward_p[t, s_ill] += log(observation_model[s_ill, o_types[id]]);
      }
      if(o_types[id] == o_pos) {
        target += normal_lpdf(viral_load[id] | viral_load_intercept + base_linpred, viral_load_sigma);
      }
    }
    target += log_sum_exp(log_forward_p[patient_max_time[p],]);

    // for(t_ in 1:(N_time - 1)) {
    //   int t_b = N_time - t_;
    //   int obs = observations[p, t_b + 1];
    //   log_backward_p[t_b, s_healthy] = log_backward_p[t_b + 1, s_healthy] + log(observation_model[s_healthy, obs]);
    //   
    //   acc_ill[s_healthy] = log_backward_p[t_b + 1, s_healthy] + log(p_ill_to_healthy) + log(observation_model[s_healthy, obs]);
    //   acc_ill[s_ill] = log_backward_p[t_b + 1, s_ill] + log1m(p_ill_to_healthy) + log(observation_model[s_ill, obs]);
    //   log_backward_p[t_b, s_ill] = log_sum_exp(acc_ill);
    // }
    // target += log_backward_p[1, s_ill]; // Initial state is always "ill"
  }

  beta ~ normal(0, fixed_prior_sd);
  time_effect ~ normal(0, time_effect_prior_sd);
  state_intercept ~ normal(0, state_intercept_prior_sd);
  viral_load_intercept ~ normal(viral_load_intercept_prior_mean, viral_load_intercept_prior_sd);
  viral_load_sigma ~ normal(0, viral_load_sigma_prior_sd);

  // for(s in 1:N_states) {
  //   target += dirichlet_lpdf(observation_model[s] | observation_model_prior[s]);
  // }
}

