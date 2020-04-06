functions {
  //TODO: Delegate to the Stan implementation
  real log_inv_logit_diff(real x, real y) {
    return(x - log1p_exp(x) + log1m_exp(y - x) - log1p_exp(y));
  }
  
  vector ordered_logistic_log_probs(real lambda, vector c) {
    int N = rows(c) + 1;
    vector[N] log_probs;
    log_probs[1] = - log1p_exp(lambda - c[1]);
    for(i in 2:(N-1)) {
      log_probs[i] = log_inv_logit_diff(lambda - c[i - 1], lambda - c[i]);
    }
    log_probs[N] = - log1p_exp(c[N - 1] - lambda);
    
    return log_probs;
  }
}

data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int<lower=1> N_time;
  int<lower=0> N_fixed;
  int<lower=1> N_ill_states;
  
  int<lower=1, upper = N_ill_states> central_ill_state;
  vector<lower=0>[N_ill_states] ill_mean_viral_load;
  
  //Constants for observation types
  int<lower=1> o_neg;
  int<lower=1> o_pos;
  int<lower=1> o_severe;
  int<lower=0, upper=1> use_severe_state;
  
  int<lower=1, upper=3> o_types[N_obs];
  int<lower=1, upper=N_patients> patients[N_obs];
  int<lower=1, upper = N_time> times[N_obs];
  vector[N_obs] viral_load;
  int<lower=0, upper=1> viral_load_known[N_obs];
  
  int N_predictor_sets;
  int<lower=1,upper=N_predictor_sets> X_index[N_patients, N_time];
  
  matrix[N_predictor_sets, N_fixed] X;
  
  
  vector<lower=0>[N_fixed] fixed_prior_sd;
  real<lower=0> observation_sigma_prior_sd;
  real<lower=0> transition_thresholds_prior_sd;
  real<lower=0> state_intercept_prior_sd;

  //ordered[N_ill_states + use_severe_state] transition_thresholds;

}

transformed data {
  int N_states = N_ill_states + 1 + use_severe_state; //ill + healthy + severe
  int s_healthy = 1;
  int ill_states_shift = 1;
  int ill_states[N_ill_states];
  int s_severe = use_severe_state ? N_states : -1;

  //0 means no observation
  int<lower=0,upper=N_obs> obs_ids[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0> patient_max_time[N_patients];
  int<lower=0> patient_min_time_severe[N_patients] = rep_array(N_time + 1, N_patients);

  for(n in 1:N_obs) {
    patient_max_time[patients[n]] = max(times[n], patient_max_time[patients[n]]);
    obs_ids[patients[n], times[n]] = n;
    if(o_types[n] == o_severe) {
      if(!use_severe_state) {
        reject("Severe observation is not allowed when use_severe_state = 0");
      }
      patient_min_time_severe[patients[n]] = min(times[n], patient_min_time_severe[patients[n]]);
    }
  }

  // Remove all but the first observation of severe (it does not inform the current model in any way)
  for(p in 1:N_patients) {
    if(patient_min_time_severe[p] < N_time) {
      for(t in patient_min_time_severe[p]:N_time) {
        obs_ids[p, t] = 0;
      }
    }
  }
  

  for(i in 1:N_ill_states) {
    ill_states[i] = i + ill_states_shift;
  }
}

parameters {
  vector[N_fixed] beta;

  ordered[N_states - 1] transition_thresholds;
  //Defines the transition model for ill states wrt. linear predictor
  positive_ordered[N_ill_states - central_ill_state] state_intercepts_high;
  positive_ordered[central_ill_state - 1] neg_state_intercepts_low;

  //base for the observation model
  real<lower=0.5, upper=1> sensitivity; 
  real<lower=0.5, upper=1> specificity; 
  real<lower=0> observation_sigma;
  
}

transformed parameters {
  ordered[N_ill_states] state_intercepts;
  
  state_intercepts[central_ill_state] = 0;
  state_intercepts[(central_ill_state + 1) : N_ill_states] = state_intercepts_high;
  for(i in 1:(central_ill_state - 1)) {
    state_intercepts[i] = -neg_state_intercepts_low[central_ill_state - i];
  }
  // simplex[N_obs_types] observation_model[N_states];
  // observation_model[s_healthy, o_neg] = specificity; 
  // observation_model[s_healthy, o_pos] = 1 - specificity;
  // observation_model[s_ill, o_neg] = 1 - sensitivity; 
  // observation_model[s_ill, o_pos] = sensitivity;
}

model {
  //Assumming ill at t = 0. Time is shifted by one, to include 0 at the start
  matrix[N_time + 1, N_states] log_forward_p;
  row_vector[N_states] observation_model_negative;
  row_vector[N_states] observation_model_positive_unknown;
  
  
  //transitions are only possible _from_ ill states
  real acc_transition[N_ill_states];
  vector[N_states] transition_log_p[N_predictor_sets, N_ill_states]; 

  //Precompute all the transition matrices
  for(x_id in 1:N_predictor_sets) {
    real linpred = 0;
    if(N_fixed > 0) {
      linpred += X[x_id, ] * beta;
    }
    for(s_index in 1:N_ill_states) {
      transition_log_p[x_id, s_index, ] = ordered_logistic_log_probs(state_intercepts[s_index] + linpred, transition_thresholds);
    }
  }


  // Init HMM, at t == 0 the state is equally separated among all ill states
  log_forward_p[1, ill_states] = rep_row_vector(-log(N_ill_states), N_ill_states);
  log_forward_p[1, s_healthy] = log(0);
  if(use_severe_state) {
    log_forward_p[1, s_severe] = log(0);
  }


  //Init constant parts of the observation model
  observation_model_negative[s_healthy] = log(specificity);
  if(use_severe_state) {
    observation_model_negative[s_severe] = log(0);
  }
  observation_model_negative[ill_states] = rep_row_vector(log1m(sensitivity), N_ill_states); //TODO consider gradual change among states

  observation_model_positive_unknown[s_healthy] =  log1m(specificity);
  if(use_severe_state) {
    observation_model_positive_unknown[s_severe] = log(0);
  }
  observation_model_positive_unknown[ill_states] = rep_row_vector(log(sensitivity), N_ill_states);

  for(p in 1:N_patients) {
    for(t in 1:patient_max_time[p]) {
      //HMM update (the HMM time is shifted by one to include t = 0)
      int t_hmm = t + 1;
      int id = obs_ids[p, t];
      
      //Transitions - only allowed from ill states
      for(s_to in 1:N_states) {
        for(s_from in ill_states) {
          acc_transition[s_from - ill_states_shift] = log_forward_p[t_hmm - 1, s_from] + transition_log_p[X_index[p,t], s_from - ill_states_shift, s_to];
        }
        if(s_to == s_healthy || (use_severe_state && s_to == s_severe)) {
          //Add the probability of already being in one of the "terminal" states
          log_forward_p[t_hmm, s_to] = log_sum_exp(append_array(acc_transition, {log_forward_p[t_hmm - 1, s_to]}));
        } else {
          log_forward_p[t_hmm, s_to] = log_sum_exp(acc_transition);
        }
      }
      
      //Observations
      if(id != 0) {
        if(o_types[id] == o_severe) {
          log_forward_p[t_hmm, s_healthy] = log(0);
          for(s in ill_states) {
            log_forward_p[t_hmm, s] = log(0);
          }
        } else if(o_types[id] == o_neg) {
          log_forward_p[t_hmm, ] += observation_model_negative;
        } else if(o_types[id] == o_pos) {
          if(viral_load_known[id]) {
            vector[N_ill_states] per_state_lp;
            if(use_severe_state) {
              log_forward_p[t_hmm, s_severe] = log(0);
            }
            for(s_index in 1:N_ill_states) {
              per_state_lp[s_index] = normal_lpdf(viral_load[id] | ill_mean_viral_load[s_index], observation_sigma);
            }
            // Healthy -> any positive state is equally probable
            log_forward_p[t_hmm, s_healthy] += log1m(specificity) -log(N_ill_states); 
            if(N_ill_states > 1) {
              log_forward_p[t_hmm, s_healthy] += log_sum_exp(per_state_lp);
            } else {
              log_forward_p[t_hmm, s_healthy] += per_state_lp[1];
            }
              
            for(s in ill_states) {
              log_forward_p[t_hmm, s] += log(sensitivity) + per_state_lp[s - ill_states_shift];
            }
          } else {
            log_forward_p[t_hmm, ] += observation_model_positive_unknown;
          }
        } else {
          reject("Unknown observation type");
        }
        
      }
    }
    target += log_sum_exp(log_forward_p[patient_max_time[p] + 1,]);
  }

  beta ~ normal(0, fixed_prior_sd);
  observation_sigma ~ normal(0, observation_sigma_prior_sd);
  transition_thresholds ~ normal(0, transition_thresholds_prior_sd);
  state_intercepts_high ~ normal(0, state_intercept_prior_sd);
  neg_state_intercepts_low ~ normal(0, state_intercept_prior_sd);
}

