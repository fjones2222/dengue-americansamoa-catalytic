//American Samoa serosurvey
//Residence history in AS, endemic, or non-endemic countries used to determine hazard

data {
  int <lower=0> N; //the number of individuals in the serosurvey
  int <lower=0> A_sp; //number of age classes in the serosurvey
  int <lower=0,upper=1> Y[N]; //serostatus of each individual
  real<lower=0> w[N]; //survey weight for each individual
  real <lower=-1,upper=2> endres_history[N,A_sp]; //the residence history of all children
  int <lower=0> r_v; //the number of covariates for lambda
  real reg_vars[N,r_v]; //covariates on lambda
  int birthyear[N]; // calendar year of birth
  int <lower=1,upper=12> birthmonth[N];
  int serosurv_year; // calendar year in which the serosurvey was conducted
  int<lower=0> N_pos_control; //number of positive controls in the validation data
  int<lower=0,upper=N_pos_control> control_tp; // number of true positive tests in the validation data
  int<lower=0> N_neg_control; // number of negative controls in the validation data
  int<lower=0,upper=N_neg_control> control_fp;// number of false positives by the diagnostic test in the validation study
  int<lower=2,upper=4>n_serotypes; //how many "endemic" serotypes circulate
  
  int <lower=0> AG; //the number of age classes
  int <lower=0> max_A; //the max age
  int casedata_lastyear; //the latest year for which case data is available
  int casedata_firstyear; // the earliest year for which case data is available
  int <lower=0>	secondary_cases[AG, casedata_lastyear-casedata_firstyear+1]; //the number of cases reported at each time point for each age group
  int <lower=0> pop[AG, casedata_lastyear-casedata_firstyear+1]; // the total population at each time point for each age group
  int <lower=0> lr_bound[AG]; // lower age groups bounds
  int <lower=0> ur_bound[AG]; // upper age groups bounds
  real foi_mean;// hyperparameter lambda mean prior
  real foi_sd;// hyperparameter lambda sd prior
  real reporting_mean;// hyperparameter reporting rate mean prior
  real reporting_sd;// hyperparameter reporting rate sd prior
  
  int old_serosurv_year; // the year of the previous serosurvey
  int old_serosurv_maxage; // the oldest age included in the previous serosurvey. 
  // Note that if old_serosurv_year - old_serosurv_maxage < casedata_firstyear - max_A, we won't have reconstructed enough
  // historical transmission to inform this seroprevalence. In this case, need to rework the code to have lambda go back to
  // min(old_serosurv_year - old_serosurv_maxage,casedata_firstyear - max_A)
  
  int pos2010_18to25;
  int total2010_18to25;
  int pos2010_26to40;
  int total2010_26to40;
  
  int clust[N]; //cluster of each individual (numeric from 1 to num_clust, must be ordered)
  int <lower=0> num_clust; //number of clusters
  real<lower=0> alpha_priormean;
  real<lower=0> alpha_priorvar;
  
  
}

transformed data {
  
  int T_case = casedata_lastyear-casedata_firstyear+1; // number of years for which we have case data
  
  int serosurv_firstyear = min(birthyear); // earliest birth year of individuals enrolled in the serosurvey
  int serosurv_lastyear = max(birthyear); // latest birth year of individuals enrolled in the serosurvey
  int T_ss = serosurv_lastyear - serosurv_firstyear + 1; // number of years from which we are estimating lambda from the serosurvey
  
  int lastyear_annualfoi = max(serosurv_lastyear,casedata_lastyear);
  
  int firstyear_foi = min(old_serosurv_year - old_serosurv_maxage,casedata_firstyear - max_A);
  
  int T_lambda = lastyear_annualfoi - firstyear_foi + 1; // the number of years for which we will estimate
  
}

parameters {
  
  real <lower=-10,upper=-0.7> logit_lambda_end; // yearly FOI in endemic countries
  real <lower=-10,upper=-0.7> logit_lambda_nonend; //yearly FOI in non-endemic countries
  real <lower=0.5, upper=1> spec; // specificity of the diagnostic test.
  real <lower=0.5, upper=1> sens; // sensitivity of the diagnostic test.
  real <lower=-4,upper=4> beta_covs[r_v]; // beta coefficients
  real <lower=-10,upper=-0.7> log_post_lambda; // average annual lambda from the last year for which we have annual data (either from
  // serosurvey or case data), as we cannot distinguish annual FOI within those years

  real <lower=0> invroot_phi; // overdispersion of the case reporting
  real reporting_rate0; // reporting rate hyperparameter
  real reporting_rate_logit[T_case]; // time varying reporting rate
  real lambda0_logit; // average FOI over time
  real lambdaRE_logit[T_lambda]; // yearly FOI random effect (estimated from case data alone)
  real <lower=0, upper=1> alpha; // proportion of first age group experiencing lambda[t-1]
  
  real <lower=0,upper=1> reprate_pvss; // relative reporting rate for primary vs secondary infections
  
  real<lower=0,upper=3> frailty_alpha; //gamma distribution of hazard by cluster, alpha
  real<lower=0> frailty[num_clust];
  
  
}

transformed parameters {
  
  // lambda, either from case or serosurvey
  
  real <lower=0, upper=1> lambda[T_lambda]; // yearly FOI
  real <lower=0> phi;
  real reporting_rate[T_case]; // reporting rate

  phi = 1/(invroot_phi)^2;

  for (t in 1:(T_case)) {
    reporting_rate[t] = inv_logit(reporting_rate_logit[t]);
  }

  for (t in 1:(T_lambda)) {
    lambda[t] = inv_logit(lambdaRE_logit[t]);
  }

  // the rest of this code is focused around constructing individual-level cumulative FOI based on each 
  // individual's age and residence history
  real lambda_i[N]; // Lambda for each individual included in the serosurvey
  real exp_sp[N]; // Probability of seropositivity for each individual included in the serosurvey
  
  real sumw;
  real lambda_rep[T_ss]; // Population-representative FOI for AS for time period of serosurvey
  
  real lambda_temp[T_lambda]; // Temporary to store lambda/lambda_rep

  real cum_lambda[max_A, T_case]; // Cumulative lambda experienced by residents of AS
  real mono[max_A, T_case]; // Proportion monotypic
  real susc[max_A, T_case]; // Proportion susceptible
  real exp_inc[max_A, T_case]; // Expected incidence
  real exp_reported_cases[AG, T_case]; // Expected reported cases
  real exp_inc_grouped[AG, T_case]; // Expected incidence in age groups


  // build up the cumulative hazard for each individual based on their residence history
  // multiply by hr for the cluster and by hrs for covariates
  
  //Lambda by individual
  for (i in 1:N) {
    
    lambda_i[i] = 0;
      // Go year by year from the calendar year of birth to the calendar year of the serosurvey
      for (y in (birthyear[i]):serosurv_year) {
        
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
          // To pick out the right entry of lambda, which goes from the casedata_firstyear - max_A to casedata_lastyear
          lambda_i[i] = lambda_i[i] + lambda[T_lambda  - (casedata_lastyear - y)];
        }
        
        if (y>lastyear_annualfoi) {
          lambda_i[i] = lambda_i[i] + exp(log_post_lambda);
        }
      }
    
    // Expected probability of seropositivity given individual-level FOI (multiplied by HRs for covariates)
    exp_sp[i] = 1-exp(-n_serotypes * lambda_i[i] * exp(dot_product(beta_covs,reg_vars[i,]))* frailty[clust[i]]) ; 
  }

  // Now get hazard experienced by age cohorts for case data
  // First get population average lambda from serosurvey
  
  // Sum of survey weights
  sumw=0;
  for (i in 1:N) {
    sumw+=w[i];
  }

  // Get lambda representative of AS population during serosurvey period by taking inverse-weighted sum
  // of FOIs from serosurvey participants
  for (y in serosurv_firstyear:serosurv_lastyear) {

    lambda_rep[y - serosurv_firstyear + 1]=0;

    for (i in 1:N) {
      if (birthyear[i]<=y) {
        // lambda_rep is representative FOI for AS based on serosurvey participants
        lambda_rep[y - serosurv_firstyear + 1] += lambda[T_lambda - (casedata_lastyear-y)] * exp(dot_product(beta_covs,reg_vars[i,])) * frailty[clust[i]]* w[i] / sumw;
      }
    }

  }

  // For the years not covered by the serosurvey, use lambda (from case data)
  for (y in (casedata_lastyear-T_lambda+1):(serosurv_firstyear-1)) {
    lambda_temp[T_lambda - (casedata_lastyear-y)] = lambda[T_lambda - (casedata_lastyear-y)];
  }

  // For the years covered by the serosurvey, use lambda_rep (from serosurvey participants)
  for (y in serosurv_firstyear:serosurv_lastyear) {
    lambda_temp[T_lambda - (casedata_lastyear-y)] = lambda_rep[y - serosurv_firstyear + 1];
  }

  for (y in (serosurv_lastyear+1):lastyear_annualfoi) {
    lambda_temp[T_lambda - (casedata_lastyear-y)] = lambda[T_lambda - (casedata_lastyear-y)];
  }

  // Cumulative FOI by age in years (a) and year considered (t)
  cum_lambda[1, 1] = alpha * lambda_temp[T_lambda - (lastyear_annualfoi - casedata_firstyear)];
  for (a in 2:max_A) {
    cum_lambda[a, 1] = sum(lambda_temp[(T_lambda - (lastyear_annualfoi - casedata_firstyear) - a + 1):(T_lambda - (lastyear_annualfoi - casedata_firstyear))]);
  }
  for (t in 2:T_case) {
    cum_lambda[1, t] = alpha * lambda_temp[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1];//first age group born with no past exposure
    for (a in 2:max_A) {
      cum_lambda[a, t] = cum_lambda[a-1, t-1] + lambda_temp[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1];
    }
  }

  // Expected case numbers (representing secondary cases)
  for (t in 1:T_case) {
    for (a in 1:max_A) {
      susc[a, t] = exp(-n_serotypes * cum_lambda[a, t]) ;
      mono[a, t] = n_serotypes * exp(-(n_serotypes-1) * cum_lambda[a, t]) * (1 - exp(-cum_lambda[a, t])) ;
      exp_inc[a, t] = reprate_pvss * n_serotypes * lambda_temp[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1] * susc[a,t] + (1-reprate_pvss) * (n_serotypes-1) * lambda_temp[T_lambda - (lastyear_annualfoi - casedata_firstyear)+t-1] * mono[a, t];
    }
  }

  // Expected reported case numbers by age groups in data (as expectation of negative binomial for the likelihood of observed case numbers)
  for (t in 1:T_case) {
    for (a_gr in 1:AG) {
      exp_inc_grouped[a_gr, t] = mean(exp_inc[lr_bound[a_gr]:ur_bound[a_gr], t]);
      exp_reported_cases[a_gr, t] = exp_inc_grouped[a_gr, t] * reporting_rate[t] * pop[a_gr, t];
    }
  }

  // Now build up the seroprevalence for the old serosurvey
  real cum_lambda_oldserosurv[old_serosurv_maxage]; // Cumulative lambda experienced by residents of AS
  real susc_oldserosurv[old_serosurv_maxage]; // Proportion susceptible

  // Cumulative FOI by age in years (a) and year considered (t)
  cum_lambda_oldserosurv[1] = alpha * lambda_temp[T_lambda - (lastyear_annualfoi - old_serosurv_year)];
  for (a in 2:old_serosurv_maxage) {
    cum_lambda_oldserosurv[a] = sum(lambda_temp[(T_lambda - (lastyear_annualfoi - old_serosurv_year)-a+1):(T_lambda - (lastyear_annualfoi - old_serosurv_year))]);
  }

  // Expected case numbers (representing secondary cases)
  for (a in 1:old_serosurv_maxage) {
    susc_oldserosurv[a] = exp(-n_serotypes * cum_lambda_oldserosurv[a]) ;
  }


  //new seroprevalence measures from 2010 study
  real <lower=0, upper=1> exp_susc_18to25;
  real <lower=0, upper=1> exp_susc_26to40;

  exp_susc_18to25 = mean(susc_oldserosurv[18:25]);
  exp_susc_26to40 = mean(susc_oldserosurv[26:40]);
  
}

model {
  
  // Priors
  alpha ~ beta(2, 1);

  reporting_rate0 ~ normal(reporting_mean, sqrt(reporting_sd^2 / 2));
  for (t in 1:T_case) {
    reporting_rate_logit[t] ~ normal(reporting_rate0, sqrt(reporting_sd^2 / 2));
  }

  lambda0_logit ~ normal(foi_mean, sqrt(foi_sd ^2 / 2));
  // logitnorm::twCoefLogitnorm(0.05, 0.3, perc = 0.975)
  for (i in 1:T_lambda) {
    lambdaRE_logit[i] ~ normal(lambda0_logit, sqrt(foi_sd ^2 / 2));
  }
  //plot(density(boot::inv.logit(rnorm(1000, rnorm(1000,-5, 1),1))))

  invroot_phi ~ normal(0, 100);
  
  //likelihood for validation data
  target+= binomial_lpmf(control_tp | N_pos_control, sens);
  target+= binomial_lpmf(control_fp | N_neg_control, 1-spec);
  
  for (i in 1:N) {
    //likelihood for seroprevalence data
    target += w[i] * bernoulli_lpmf(Y[i] | exp_sp[i]*sens+(1-exp_sp[i])*(1-spec));
  }
  
  // likelihood for case data
  for (t in 1:T_case) {
    for (a_gr in 1:AG) {
      secondary_cases[a_gr, t] ~ neg_binomial_2(exp_reported_cases[a_gr, t], phi);
    }
  }
  
  // likelihood for new seroprevalence data

  target+= binomial_lpmf(total2010_18to25-pos2010_18to25 | total2010_18to25, exp_susc_18to25);
  target+= binomial_lpmf(total2010_26to40-pos2010_26to40| total2010_26to40, exp_susc_26to40);

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
