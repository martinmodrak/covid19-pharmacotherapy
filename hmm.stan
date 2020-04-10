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
  
  //HMM functions
  row_vector init_observation_model_negative(
    int N_states, int[] ill_states, int s_healthy, int s_severe, int use_severe_state, 
    real log_specificity, real log1m_sensitivity
  ) {
    row_vector[N_states] observation_model_negative;
    observation_model_negative[s_healthy] = log_specificity;
    if(use_severe_state) {
      observation_model_negative[s_severe] = log(0);
    }
    observation_model_negative[ill_states] = rep_row_vector(log1m_sensitivity, size(ill_states)); //TODO consider gradual change among states
    
    return observation_model_negative;
  }

  row_vector init_observation_model_positive_unknown(
    int N_states, int[] ill_states, int s_healthy, int s_severe, int use_severe_state, 
    real log1m_specificity, real log_sensitivity
  ) {
    row_vector[N_states] observation_model_positive_unknown;
    observation_model_positive_unknown[s_healthy] =  log1m_specificity;
    if(use_severe_state) {
      observation_model_positive_unknown[s_severe] = log(0);
    }
    observation_model_positive_unknown[ill_states] = rep_row_vector(log_sensitivity, size(ill_states));
    
    return observation_model_positive_unknown;
  }
  
  // Compute a single transition matrix
  // The indices are [to_state, from_state] 
  // where to_state range 1 .. N_states and from_state 1 .. N_ill_states
  matrix compute_transition_log_p(
    int N_states, int N_ill_states, int N_fixed, row_vector X,
    vector beta, vector transition_thresholds, vector state_intercepts
  ) {
      matrix[N_states, N_ill_states] transition_log_p;
      real linpred = 0;
      if(N_fixed > 0) {
        linpred += X * beta;
      }
      for(s_index in 1:N_ill_states) {
        transition_log_p[, s_index] = ordered_logistic_log_probs(state_intercepts[s_index] + linpred, transition_thresholds);
      }
      return transition_log_p;
  }

  //Precompute all the transition matrices
  //The indices are [predictor_set, to_state, from_state] 
  //where to_state range 1 .. N_states and from_state 1 .. N_ill_states
  matrix[] compute_all_transition_log_p(
    int N_states, int N_predictor_sets, int N_ill_states, int N_fixed, matrix X,
    vector beta, vector transition_thresholds, vector state_intercepts
  ) {
    //transitions are only possible _from_ ill states
    matrix[N_states, N_ill_states] all_transition_log_p[N_predictor_sets]; 
  
    for(x_id in 1:N_predictor_sets) {
      all_transition_log_p[x_id] = compute_transition_log_p(N_states, N_ill_states, N_fixed,
        X[x_id,], beta, transition_thresholds, state_intercepts);
    }
    return all_transition_log_p;
  }
  
  matrix forward_pass(
    int patient_id, 
    int N_time, int observation_shift, int N_states, 
    int[] ill_states, int ill_states_shift, 
    vector ill_mean_viral_load, 
    int s_healthy, int s_severe, int use_severe_state, 
    int[] X_index_patient, int[] o_types, 
    int[,] obs_ids,
    int o_neg, int o_pos, int o_severe,
    int[] viral_load_known, vector viral_load,
    real log1m_specificity, real log_sensitivity, real observation_sigma,
    matrix[] all_transition_log_p, 
    row_vector observation_model_negative, row_vector observation_model_positive_unknown
  ) {
    int N_ill_states = size(ill_states);
    //Assumming ill at t = 0. Time is shifted by one, to include 0 at the start
    matrix[N_time + 1, N_states] log_forward_p;
    //transitions are only possible _from_ ill states
    real acc_transition[N_ill_states];
    
    
    // Init HMM, at t == 0 the state is equally likely to be any of the ill states
    log_forward_p[1, ill_states] = rep_row_vector(-log(N_ill_states), N_ill_states);
    log_forward_p[1, s_healthy] = log(0);
    if(use_severe_state) {
      log_forward_p[1, s_severe] = log(0);
    }
    
    for(t in 1:N_time) {
      //HMM update (the HMM time is shifted by one to include t = 0)
      int t_hmm = t + 1;
      int o_id;
      if(t - observation_shift > 0) {
        o_id = obs_ids[patient_id, t - observation_shift];
      } else {
        o_id = 0;
      }

      //Transitions - only allowed from ill states
      for(s_to in 1:N_states) {
        for(s_from in ill_states) {
          acc_transition[s_from - ill_states_shift] = log_forward_p[t_hmm - 1, s_from] + all_transition_log_p[X_index_patient[t], s_to, s_from - ill_states_shift];
        }
        if(s_to == s_healthy || (use_severe_state && s_to == s_severe)) {
          //Add the probability of already being in one of the "terminal" states
          log_forward_p[t_hmm, s_to] = log_sum_exp(append_array(acc_transition, {log_forward_p[t_hmm - 1, s_to]}));
        } else {
          log_forward_p[t_hmm, s_to] = log_sum_exp(acc_transition);
        }
      }
      
      //Observations
      if(o_id != 0) {
        if(o_types[o_id] == o_severe) {
          log_forward_p[t_hmm, s_healthy] = log(0);
          for(s in ill_states) {
            log_forward_p[t_hmm, s] = log(0);
          }
        } else if(o_types[o_id] == o_neg) {
          log_forward_p[t_hmm, ] += observation_model_negative;
        } else if(o_types[o_id] == o_pos) {
          if(viral_load_known[o_id] && N_ill_states > 1) {
            vector[N_ill_states] per_state_lp;
            if(use_severe_state) {
              log_forward_p[t_hmm, s_severe] = log(0);
            }
            for(s_index in 1:N_ill_states) {
              //Normal truncated at 0
              per_state_lp[s_index] = 
                normal_lpdf(viral_load[o_id] | ill_mean_viral_load[s_index], observation_sigma)
                -normal_lccdf(0 | ill_mean_viral_load[s_index], observation_sigma);
            }
            // Healthy -> any positive state is equally probable
            log_forward_p[t_hmm, s_healthy] += log1m_specificity -log(N_ill_states); 
            if(N_ill_states > 1) {
              log_forward_p[t_hmm, s_healthy] += log_sum_exp(per_state_lp);
            } else {
              log_forward_p[t_hmm, s_healthy] += per_state_lp[1];
            }
              
            for(s in ill_states) {
              log_forward_p[t_hmm, s] += log_sensitivity + per_state_lp[s - ill_states_shift];
            }
          } else {
            log_forward_p[t_hmm, ] += observation_model_positive_unknown;
          }
        } else {
          reject("Unknown observation type");
        }
        
      }
      
    }
    
    return log_forward_p; 
  }
}

data {
  int<lower=1> N_patients;
  int<lower=1> N_obs;
  int<lower=1> N_time;
  int<lower=0> N_fixed;
  int<lower=1> N_ill_states;
  
  //Central state will have intercept = 0
  int<lower=1, upper = N_ill_states> central_ill_state;
  vector<lower=0>[N_ill_states] ill_mean_viral_load;
  
  // Patients with unknown start time
  int<lower=0> N_unknown_shift;
  int<lower=1, upper = N_patients> unknown_shift_patients[N_unknown_shift];
  int<lower=0> max_observation_shift;
  
  
  // Constants for observation types, passed along to make sure
  // both R and Stan use them consistently
  int<lower=1> o_neg;
  int<lower=1> o_pos;
  int<lower=1> o_severe;
  int<lower=0, upper=1> use_severe_state; // Switches of the severe state if necessary
  
  // Observation data in long format
  int<lower=1, upper=3> o_types[N_obs];
  int<lower=1, upper=N_patients> patients[N_obs];
  int<lower=1, upper = N_time> times[N_obs];
  vector[N_obs] viral_load;
  // Whether viral load is known 
  int<lower=0, upper=1> viral_load_known[N_obs];
  
  // Note: predictors are needed even for the times when nothing is observed
  //
  // Trick for efficiency - since the predictors don't vary that much
  // I provide indices to a smaller X array. This transformation is done in R using the 
  // transform_predictors_to_unique function as it is easier there
  int N_predictor_sets;
  int<lower=1,upper=N_predictor_sets> X_index[N_patients, N_time];
  // Dirty trick: since time can potentially be an effect, I let the user specify complete
  // predictor sequences for each time point and potential shift due to unknown symptom onset,
  // As it is easier to just do the work in R.
  // (once again indexing the small X array)
  int<lower=1,upper=N_predictor_sets> X_unknown_shift_index[N_unknown_shift, N_unknown_shift > 0 ? max_observation_shift + 1 : 0, N_time];
  
  // The unique predictor sets among all patient x time points.
  matrix[N_predictor_sets, N_fixed] X;
  
  // Prior parameters
  vector<lower=0>[N_fixed] fixed_prior_sd;
  real<lower=0> observation_sigma_prior_sd;
  real<lower=0> transition_thresholds_prior_sd;
  real<lower=0> state_intercept_prior_sd;

  // Parameters for prediction
  int<lower=0, upper=1> generate_predictions;
  int<lower=0> N_new_predictions;
  matrix[N_time, N_fixed] X_new_prediction[N_new_predictions];
}

transformed data {
  int N_states = N_ill_states + 1 + use_severe_state; //ill + healthy + severe
  
  // Constants to be able to adress states by name
  int s_healthy = 1;
  int ill_states_shift = 1;
  int ill_states[N_ill_states];
  int s_severe = use_severe_state ? N_states : -1;

  //Map from patient and time to the ragged observation data
  // 0 means no observation
  int<lower=0,upper=N_obs> obs_ids[N_patients, N_time] = rep_array(0, N_patients, N_time);
  int<lower=0> patient_max_time[N_patients];
  int<lower=0> patient_min_time_severe[N_patients] = rep_array(N_time + 1, N_patients);
  // Map from patients to arrays with information about unknown symptom onset
  // 0 means known start
  int<lower=0,upper = N_unknown_shift> unknown_shift_ids[N_patients] = rep_array(0, N_patients);

  int N_time_for_prediction = (generate_predictions ? N_time : 0);

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

  // Remove all but the first observation of severe 
  // It does not inform the current model in any way, but breaks gradients
  for(p in 1:N_patients) {
    if(patient_min_time_severe[p] < N_time) {
      for(t in (patient_min_time_severe[p] + 1):N_time) {
        obs_ids[p, t] = 0;
      }
      patient_max_time[p] = patient_min_time_severe[p];
    }
  }
  
  for(up in 1:N_unknown_shift) {
    unknown_shift_ids[unknown_shift_patients[up]] = up;
    if(patient_max_time[unknown_shift_patients[up]] + max_observation_shift > N_time) {
      reject("N_time must be large enough to incorporate latest observation from each patient with unknown start and the time shift");
    }
  }

  for(i in 1:N_ill_states) {
    ill_states[i] = i + ill_states_shift;
  }
}

parameters {
  vector[N_fixed] beta;

  // Defines the transition model for ill states wrt. linear predictor
  // Threshold for the cumulative logit
  ordered[N_states - 1] transition_thresholds; 
  // State intercepts above the central state (which has 0 as intercept)
  positive_ordered[N_ill_states - central_ill_state] state_intercepts_high;
  // State intercepts below the central state
  positive_ordered[central_ill_state - 1] neg_state_intercepts_low;

  // Base for the observation model
  real<lower=0.5, upper=1> sensitivity; 
  real<lower=0.5, upper=1> specificity; 
  real<lower=0> observation_sigma;
 
  // Marginalizing the unknown start times
  simplex[max_observation_shift + 1] obs_shift_probs[N_unknown_shift]; 
}

transformed parameters {
  ordered[N_ill_states] state_intercepts;
  
  state_intercepts[central_ill_state] = 0;
  state_intercepts[(central_ill_state + 1) : N_ill_states] = state_intercepts_high;
  for(i in 1:(central_ill_state - 1)) {
    state_intercepts[i] = -neg_state_intercepts_low[central_ill_state - i];
  }
}

model {
  real log_specificity = log(specificity);
  real log1m_specificity = log1m(specificity);
  real log_sensitivity = log(sensitivity);
  real log1m_sensitivity = log1m(sensitivity);
  
  
  //Init constant parts of the observation model
  row_vector[N_states] observation_model_negative = 
    init_observation_model_negative(N_states, ill_states, s_healthy, s_severe, use_severe_state, 
    log_specificity, log1m_sensitivity);
  row_vector[N_states] observation_model_positive_unknown = 
    init_observation_model_positive_unknown(N_states, ill_states, s_healthy, s_severe, use_severe_state, 
    log1m_specificity, log_sensitivity);
  
  
  //Precompute all the transition matrices for all unique predictor sets
  matrix[N_states, N_ill_states] all_transition_log_p[N_predictor_sets] = 
    compute_all_transition_log_p(
      N_states, N_predictor_sets, N_ill_states, N_fixed, X,
      beta, transition_thresholds, state_intercepts
    );

  for(p in 1:N_patients) {
    if(unknown_shift_ids[p] == 0) {
      //Known start time
      int t_max = patient_max_time[p];
      int observation_shift = 0;
      int X_index_patient[N_time] = X_index[p,];
      matrix[t_max + 1, N_states] log_forward_p =
        forward_pass(
          p, 
          t_max, observation_shift, N_states, 
          ill_states, ill_states_shift, 
          ill_mean_viral_load, 
          s_healthy, s_severe, use_severe_state, 
          X_index_patient, o_types, 
          obs_ids,
          o_neg, o_pos, o_severe,
          viral_load_known, viral_load,
          log1m_specificity, log_sensitivity, observation_sigma,
          all_transition_log_p, 
          observation_model_negative, observation_model_positive_unknown
      );
      
      target += log_sum_exp(log_forward_p[t_max + 1,]);
    } else {
      // Unknown start time - marginalize over all possible start times
      vector[max_observation_shift + 1] per_shift_lp;
      for(observation_shift in 0:max_observation_shift) {
        int t_max = patient_max_time[p] + observation_shift;
        int X_index_patient[N_time] = X_unknown_shift_index[unknown_shift_ids[p], observation_shift + 1, ];
        matrix[t_max + 1, N_states] log_forward_p =
          forward_pass(
            p, 
            t_max, observation_shift, N_states, 
            ill_states, ill_states_shift, 
            ill_mean_viral_load, 
            s_healthy, s_severe, use_severe_state, 
            X_index_patient, o_types, 
            obs_ids,
            o_neg, o_pos, o_severe,
            viral_load_known, viral_load,
            log1m_specificity, log_sensitivity, observation_sigma,
            all_transition_log_p, 
            observation_model_negative, observation_model_positive_unknown
        );
        
        per_shift_lp[observation_shift + 1] = 
          log(obs_shift_probs[unknown_shift_ids[p], observation_shift + 1]) + log_sum_exp(log_forward_p[t_max + 1,]);
        
      }
      
      target += log_sum_exp(per_shift_lp);
    }
  }

  beta ~ normal(0, fixed_prior_sd);
  observation_sigma ~ normal(0, observation_sigma_prior_sd);
  transition_thresholds ~ normal(0, transition_thresholds_prior_sd);
  state_intercepts_high ~ normal(0, state_intercept_prior_sd);
  neg_state_intercepts_low ~ normal(0, state_intercept_prior_sd);
}

generated quantities {
  int<lower=0, upper=max_observation_shift> observation_shift_pred[N_unknown_shift];
  int<lower=1, upper=N_states> state_pred_new[N_new_predictions, N_time_for_prediction];
  real log_p_new[N_new_predictions, N_time_for_prediction, N_states];
  int<lower=1, upper=N_states> state_pred[N_patients, N_time_for_prediction];
  int<lower=1> o_type_pred[N_patients, N_time_for_prediction];
  matrix[N_patients, N_time_for_prediction] viral_load_pred;
  
  for(up in 1:N_unknown_shift) {
    observation_shift_pred[up] = categorical_rng(obs_shift_probs[up]) - 1;
  }

  if(generate_predictions) {
    real log_specificity = log(specificity);
    real log1m_specificity = log1m(specificity);
    real log_sensitivity = log(sensitivity);
    real log1m_sensitivity = log1m(sensitivity);
    
    //Init constant parts of the observation model
    row_vector[N_states] observation_model_negative = 
      init_observation_model_negative(N_states, ill_states, s_healthy, s_severe, use_severe_state, 
      log_specificity, log1m_sensitivity);
    row_vector[N_states] observation_model_positive_unknown = 
      init_observation_model_positive_unknown(N_states, ill_states, s_healthy, s_severe, use_severe_state, 
      log1m_specificity, log_sensitivity);
    
    
    //Precompute all the transition matrices
    matrix[N_states, N_ill_states] all_transition_log_p[N_predictor_sets] = 
      compute_all_transition_log_p(
        N_states, N_predictor_sets, N_ill_states, N_fixed, X,
        beta, transition_thresholds, state_intercepts
      );
      
    for(p in 1:N_patients) {
      int observation_shift;
      int X_index_patient[N_time];
      
      if(unknown_shift_ids[p] == 0) {
        observation_shift = 0;
        X_index_patient = X_index[p, ];
      } else {
        observation_shift = observation_shift_pred[unknown_shift_ids[p]];
        X_index_patient = X_unknown_shift_index[unknown_shift_ids[p], observation_shift + 1, ];
      }
      {
          matrix[N_time_for_prediction + 1, N_states] log_forward_p;

          log_forward_p =
            forward_pass(
              p, 
              N_time_for_prediction, observation_shift, N_states, 
              ill_states, ill_states_shift, 
              ill_mean_viral_load, 
              s_healthy, s_severe, use_severe_state, 
              X_index_patient, o_types, 
              obs_ids,
              o_neg, o_pos, o_severe,
              viral_load_known, viral_load,
              log1m_specificity, log_sensitivity, observation_sigma,
              all_transition_log_p, 
              observation_model_negative, observation_model_positive_unknown
          );
          
          //Backward sampling pass
          
          state_pred[p, N_time_for_prediction] = 
            categorical_rng(softmax(to_vector(log_forward_p[N_time_for_prediction + 1,])));
          //Backward pass to generate samples
          for(t_inv in 2:N_time_for_prediction) {
            int t = N_time_for_prediction - t_inv + 1;
            int t_hmm = t + 1;
            int next_state = state_pred[p, t + 1];
            
            vector[N_states] log_probs;
            
            log_probs[ill_states] = to_vector(
              log_forward_p[t_hmm, ill_states] + all_transition_log_p[X_index[p, t], next_state,]);
              
            if(next_state == s_healthy) {
              log_probs[s_healthy] = log_forward_p[t_hmm, s_healthy];
              if(use_severe_state) {
                log_probs[s_severe] = log(0);
              }
            } else if(next_state == s_severe) {
              log_probs[s_healthy] = log(0);
              log_probs[s_severe] = log_forward_p[t_hmm, s_severe];
            } else {
              log_probs[s_healthy] = log(0);
              if(use_severe_state) {
                log_probs[s_severe] = log(0);
              }
            }
            
            state_pred[p, t] = categorical_rng(softmax(log_probs));
          }
  
          //Simulate observations from hidden states
          for(t in 1:N_time_for_prediction) {
            if(state_pred[p, t] == s_severe) {
              o_type_pred[p, t] = o_severe;
              viral_load_pred[p, t] = 0;
            } else if(state_pred[p, t] == s_healthy) {
              if(bernoulli_rng(specificity)) {
                o_type_pred[p, t] = o_neg;
                viral_load_pred[p, t] = 0;
              } else {
                real mean_load = ill_mean_viral_load[categorical_rng(rep_vector(inv(N_ill_states), N_ill_states))];
                o_type_pred[p, t] = o_pos;
                viral_load_pred[p, t] = 0;
                while(viral_load_pred[p, t] <= 0) {
                  viral_load_pred[p, t] = normal_rng(mean_load, observation_sigma);
                }
              }
            } else {
              if(bernoulli_rng(sensitivity)) {
                real mean_load = ill_mean_viral_load[state_pred[p, t] - ill_states_shift];
                o_type_pred[p, t] = o_pos;
                viral_load_pred[p, t] = 0;
                while(viral_load_pred[p, t] <= 0) {
                  viral_load_pred[p, t] = normal_rng(mean_load, observation_sigma);
                }
              } else {
                o_type_pred[p, t] = o_neg;
                viral_load_pred[p, t] = 0;
              }
            }
          }
          
        } //log_forward_p scope
      } //cycle over patients
    
      //New patient predictions (ignoring observed data)
      for(np in 1:N_new_predictions) {
        int state = categorical_rng(rep_vector(inv(N_ill_states), N_ill_states)) + ill_states_shift;
        row_vector[N_states] log_p_last;
        log_p_last[ill_states] = rep_row_vector(-log(N_ill_states), N_ill_states);
        log_p_last[s_healthy] = log(0);
        if(use_severe_state) {
          log_p_last[s_severe] = log(0);
        }
        for(t in 1:N_time) {
          real acc_transition[N_ill_states];
          matrix[N_states, N_ill_states] transition_log_p = 
            compute_transition_log_p(
              N_states, N_ill_states, N_fixed, 
              X_new_prediction[np, t, ], beta, transition_thresholds, state_intercepts);
              
          //Transitions - only allowed from ill states
          for(s_to in 1:N_states) {
            for(s_from in ill_states) {
              acc_transition[s_from - ill_states_shift] = log_p_last[s_from] + transition_log_p[s_to, s_from - ill_states_shift];
            }
            if(s_to == s_healthy || (use_severe_state && s_to == s_severe)) {
              //Add the probability of already being in one of the "terminal" states
              log_p_new[np, t, s_to] = log_sum_exp(append_array(acc_transition, {log_p_last[s_to]}));
            } else {
              log_p_new[np, t, s_to] = log_sum_exp(acc_transition);
            }
          }
          log_p_last = to_row_vector(log_p_new[np, t, ]);

          if(state != s_healthy && state != s_severe) {
            state = categorical_rng(exp(transition_log_p[, state - ill_states_shift]));
          }
          state_pred_new[np, t] = state;
        }
      }
      
  } //if generate predictions
}
