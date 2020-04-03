functions {
  // Code for Kalman filter based on the implementation in
  // https://github.com/stan-dev/math/blob/develop/stan/math/prim/prob/gaussian_dlm_obs_lpdf.hpp
  // Using the same notation, simplifed because both design matrix (F) and transition matrix (G) are
  // identity

  
  //Kalman state = [m, C, lp], 
  //m - estimated mean vector, C estimated covariance
  //lp cumulated log probability
  vector init_kalman_state(real initial_viral_load_mu, real initial_viral_load_sd) {
    return to_vector({initial_viral_load_mu, initial_viral_load_sd ^ 2, 0});
  }
  
  int kalman_state_length() {
    return 3;
  }
  
  vector kalman_step(vector kalman_state, real x, real process_variance, real obs_variance, real step_size) {
    real m = kalman_state[1];
    real C = kalman_state[2];
    real lp = kalman_state[3];
    real Q;
    real Q_inv;
    real e;
    real A;
    C += process_variance * step_size;
    //f = m;
    Q = C + obs_variance;
    Q_inv = inv(Q);
    A = C * Q_inv;
    //C -= Q * A * A;
    C *= 1 - A;
    e = x - m;
    m += A*e;
    lp += -0.5 * (log(Q) + e * e * Q_inv);
    return to_vector({m, C, lp});
  }
  
  real kalman_lp(vector kalman_state) {
    return kalman_state[3];
  }
}

data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int<lower=1> N_time;
  int<lower=0> N_fixed;
  
  //Constants for observation types
  int<lower=1> o_neg;
  int<lower=1> o_pos;
  
  //Observations: 1 = NEG, 2 = POS, 3 = ICU/DEATH
  int<lower=1, upper=2> o_types[N_obs];
  int<lower=1, upper=N_patients> patients[N_obs];
  int<lower=1, upper = N_time> times[N_obs];
  vector[N_obs] viral_load;
  int<lower=0, upper=1> viral_load_known[N_obs];
  matrix[N_time, N_fixed] X[N_patients];
  
  
  vector<lower=0>[N_fixed] fixed_prior_sd;
  real state_intercept_prior_mean;
  real<lower=0> state_intercept_prior_sd;
  real<lower=0> viral_load_intercept_prior_sd;

  real initial_viral_load_mu;
  real<lower=0> initial_viral_load_sd;
  
  real<lower=0> kalman_total_noise_prior_alpha;
  real<lower=0> kalman_total_noise_prior_beta;
  
}

transformed data {
  int N_obs_types = 2;
  int N_states = 2;
  int s_healthy = 1;
  int s_ill = 2;
  //0 means no observation
  int<lower=0,upper=N_obs> obs_ids[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0> patient_max_time[N_patients];
  int<lower=0, upper = N_time> patient_last_load[N_patients] = rep_array(0, N_patients);
  
  for(n in 1:N_obs) {
    patient_max_time[patients[n]] = max(times[n], patient_max_time[patients[n]]);
    obs_ids[patients[n], times[n]] = n;
    if(viral_load_known[n]) {
      patient_last_load[patients[n]] = max(times[n],patient_last_load[patients[n]]);
    }
  }
}

parameters {
  vector[N_fixed] beta;
  real state_intercept;
  real viral_load_intercept;
  real<lower=0> kalman_total_noise;
  real<lower=0,upper=1> kalman_process_noise_frac;
  
  real<lower=0.5, upper=1> sensitivity; 
  real<lower=0.5, upper=1> specificity; 
  

}

transformed parameters {
  real<lower=0> kalman_process_variance = (kalman_total_noise ^ 2) * kalman_process_noise_frac;
  real<lower=0> kalman_obs_variance = (kalman_total_noise ^ 2) * (1 - kalman_process_noise_frac);
  simplex[N_obs_types] observation_model[N_states];
  observation_model[s_healthy, o_neg] = specificity; 
  observation_model[s_healthy, o_pos] = 1 - specificity;
  observation_model[s_ill, o_neg] = 1 - sensitivity; 
  observation_model[s_ill, o_pos] = sensitivity;
}

model {
  //Assumming ill at t = 0. Time is shifted by one, to include 0 at the start
  matrix[N_time + 1, N_obs_types] log_forward_p;
  //matrix[N_time, N_obs_types] log_backward_p;
  real acc_healthy[2];
  real acc_ill[2];
  
  // Init HMM, at t == 0 the state is "ill", same for all patients
  log_forward_p[1, s_ill] = 0;
  log_forward_p[1, s_healthy] = log(0);
  
  //log_backward_p[N_time, ] = rep_row_vector(0, N_obs_types);
  for(p in 1:N_patients) {
    vector[kalman_state_length()] kalman_state;
    real cumulative_viral_load_linpred = 0;
    int kalman_step_size = 1; //Tracks number of time steps since last viral load observation
    if(patient_last_load[p] != 0) {
      //Init Kalman only for patients with at least one viral load observation
      kalman_state = init_kalman_state(initial_viral_load_mu, initial_viral_load_sd);
    }
    
    for(t in 1:patient_max_time[p]) {
      int id = obs_ids[p, t];
      real base_linpred = 0;
      if(N_fixed > 0) {
        base_linpred += X[p, t,] * beta;
      }

      //HMM update (the HMM time is shifted by one to includ t = 0)
      {
        int t_hmm = t + 1;
        real p_ill_to_healthy = inv_logit(state_intercept + base_linpred);
        acc_healthy[s_healthy] = log_forward_p[t_hmm - 1, s_healthy];
        acc_healthy[s_ill] = log_forward_p[t_hmm - 1, s_ill] + log(p_ill_to_healthy);
        log_forward_p[t_hmm, s_healthy] = log_sum_exp(acc_healthy);
        if(id != 0) {
          log_forward_p[t_hmm, s_healthy] += log(observation_model[s_healthy, o_types[id]]);
        }
        
        log_forward_p[t_hmm, s_ill] =  log_forward_p[t_hmm - 1, s_ill] + 
          log1m(p_ill_to_healthy);
        if(id != 0) {
          log_forward_p[t_hmm, s_ill] += log(observation_model[s_ill, o_types[id]]);
        }
      }
      
      //Kalman update - only if there are viral_loads to yet be observed
      if(t <= patient_last_load[p]) {
        real viral_load_linpred = viral_load_intercept + base_linpred;
        cumulative_viral_load_linpred += viral_load_linpred;
        if(id != 0 && viral_load_known[id]) {
          real viral_load_centered = viral_load[id] - cumulative_viral_load_linpred;
          kalman_state = kalman_step(kalman_state, viral_load_centered, kalman_process_variance, kalman_obs_variance, kalman_step_size);
          kalman_step_size = 1;
        } else {
          kalman_step_size += 1;
        }
      }
    }
    target += log_sum_exp(log_forward_p[patient_max_time[p] + 1,]);
    if(patient_last_load[p] != 0) {
      //Add accumulated lp from Kalman (if the patient has at least one viral load observation)
      target += kalman_lp(kalman_state);
    }

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
  state_intercept ~ normal(state_intercept_prior_mean, state_intercept_prior_sd);
  viral_load_intercept ~ normal(0, viral_load_intercept_prior_sd);
  kalman_total_noise ~ inv_gamma(kalman_total_noise_prior_alpha, kalman_total_noise_prior_beta);
  
  

  // for(s in 1:N_states) {
  //   target += dirichlet_lpdf(observation_model[s] | observation_model_prior[s]);
  // }
}

