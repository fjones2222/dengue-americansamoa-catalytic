//American Samoa serosurvey
//Residence history in AS, endemic, or non-endemic countries used to determine hazard

data {
  int <lower=0> N; //the number of individuals in the serosurvey
  int <lower=0> A_sp; //number of age classes in the serosurvey
  int <lower=0,upper=1> Y[N]; //serostatus of each individual
  int firstyear_exposure[N]; // calendar year of first dengue season experienced
  int serosurv_year; // calendar year in which the serosurvey was conducted
  
  real<lower=0> w[N]; //survey weight for each individual
  real <lower=-1,upper=2> endres_history[N,A_sp]; //the residence history of all children
  int <lower=0> r_v; //the number of covariates for lambda
  real reg_vars[N,r_v]; //covariates on lambda
  
  int<lower=0> N_pos_control; //number of positive controls in the validation data
  int<lower=0,upper=N_pos_control> control_tp; // number of true positive tests in the validation data
  int<lower=0> N_neg_control; // number of negative controls in the validation data
  int<lower=0,upper=N_neg_control> control_fp;// number of false positives by the diagnostic test in the validation study
  int<lower=2,upper=4>n_serotypes; //how many "endemic" serotypes circulate
  
  real foi_mean;// hyperparameter lambda mean prior
  real foi_sd_prior;// hyperparameter lambda sd prior
  real foi_sd; // fixed annual sd of foi (log-odds scale)
  
  int clust[N]; //cluster of each individual (numeric from 1 to num_clust, must be ordered)
  int <lower=0> num_clust; //number of clusters
  real<lower=0> alpha_priormean;
  real<lower=0> alpha_priorvar;
  
  
}

transformed data {
  
  int serosurv_firstyear = min(firstyear_exposure); // earliest birth year of individuals enrolled in the serosurvey
  int serosurv_lastyear = max(firstyear_exposure); // latest birth year of individuals enrolled in the serosurvey
  
  int T_lambda = serosurv_lastyear - serosurv_firstyear; // the number of years for which we will estimate annual FOI
  
}

parameters {
  
  real <lower=-10,upper=-0.7> logit_lambda_end; // yearly FOI in endemic countries
  real <lower=-10,upper=-0.7> logit_lambda_nonend; //yearly FOI in non-endemic countries
  real <lower=0.5, upper=1> spec; // specificity of the diagnostic test.
  real <lower=0.5, upper=1> sens; // sensitivity of the diagnostic test.
  real <lower=-4,upper=4> beta_covs[r_v]; // beta coefficients
  real <lower=-10,upper=-0.7> logit_post_lambda; // average annual lambda from the last year for which we have annual data (either from
  // serosurvey or case data), as we cannot distinguish annual FOI within those years
  
  real lambda0_logit; // average FOI over time
  real lambdaRE_logit[T_lambda]; // yearly FOI random effect (estimated from case data alone)
  
  //real<lower=0> foi_sd; // annual variation in FOI (log-odds scale)
  
  real<lower=0,upper=3> frailty_alpha; //gamma distribution of hazard by cluster, alpha
  real<lower=0> frailty[num_clust];
  
}

transformed parameters {
  
  // lambda, either from case or serosurvey
  
  real <lower=0, upper=1> lambda[T_lambda]; // yearly FOI
  
  real <lower=0, upper=1> post_lambda;
  
  post_lambda = inv_logit(logit_post_lambda);
  
  for (t in 1:T_lambda) {
    lambda[t] = inv_logit(lambdaRE_logit[t]);
  }
  
  // the rest of this code is focused around constructing individual-level cumulative FOI based on each 
  // individual's age and residence history
  real lambda_i[N]; // Lambda for each individual included in the serosurvey
  real exp_sp[N]; // Probability of seropositivity for each individual included in the serosurvey
  real exp_sn[N]; // Probability of seronegativity for each individual included in the serosurvey (calculated directly to avoid rounding issues)
  
  // build up the cumulative hazard for each individual based on their residence history
  // multiply by hr for the cluster and by hrs for covariates
  
  //Lambda by individual
  for (i in 1:N) {
    
    lambda_i[i] = 0;
      // Go year by year from the calendar year of birth to the calendar year of the serosurvey
      for (y in firstyear_exposure[i]:serosurv_year) {
        
        if (y>serosurv_lastyear-1) {
          lambda_i[i] = lambda_i[i] + post_lambda;
        } else {
        
          // a = age is y-serosurv_firstyear+1
          
          // Note this will only add to individual FOI if endres_history[i,a] is 0, 1, or 2
          // For years before child is born, endres_history is set to -1
  
          if (endres_history[i,y-serosurv_firstyear+1]==0) {
            lambda_i[i] = lambda_i[i] + inv_logit(logit_lambda_nonend);
          }
  
          if (endres_history[i,y-serosurv_firstyear+1]==1) {
            lambda_i[i] = lambda_i[i] + inv_logit(logit_lambda_end);
          }
  
          if (endres_history[i,y-serosurv_firstyear+1]==2) {
            // To pick out the right entry of lambda, which goes from the serosurv_firstyear + 1 to serosurv_lastyear
            lambda_i[i] = lambda_i[i] + lambda[T_lambda  - (serosurv_lastyear - y) + 1];
          }
        }
        
        
      }
    
    // Expected probability of seropositivity given individual-level FOI (multiplied by HRs for covariates)
    exp_sp[i] = 1-exp(-n_serotypes * lambda_i[i] * exp(dot_product(beta_covs,reg_vars[i,]))* frailty[clust[i]]) ; 
    exp_sn[i] = exp(-n_serotypes * lambda_i[i] * exp(dot_product(beta_covs,reg_vars[i,]))* frailty[clust[i]]) ; 
  }

}

model {
  
  // Priors
  lambda0_logit ~ normal(foi_mean, foi_sd_prior);
  for (i in 1:T_lambda) {
    lambdaRE_logit[i] ~ normal(lambda0_logit, foi_sd);
  }
  
  post_lambda ~ normal(foi_mean, foi_sd_prior);
  
  //likelihood for validation data
  target+= binomial_lpmf(control_tp | N_pos_control, sens);
  target+= binomial_lpmf(control_fp | N_neg_control, 1-spec);
  
  for (i in 1:N) {
    //likelihood for seroprevalence data
    target += w[i] * bernoulli_lpmf(Y[i] | exp_sp[i]*sens+(exp_sn[i])*(1-spec));
  }
  
  // Frailty by cluster in serosurvey
  frailty_alpha ~ lognormal(alpha_priormean,alpha_priorvar);
  for (i in 1:num_clust) {
    frailty[i] ~ gamma(frailty_alpha,frailty_alpha);
  }


}

generated quantities {
  real exp_obs_sp[N];
  for (i in 1:N) {
    
    // model predicted seropositivity (accounting for sensitivity and specificity)
    exp_obs_sp[i] = exp_sp[i]*sens+(1-exp_sp[i])*(1-spec);

  }
  
}
