# reporting parameter priors
reporting_mean <- -2.2
reporting_sd <- 0.7 #2 is what is used to be

## Code to get prior parameters for the yearly FOI, based on observed incidence
incidence <- calc_incidence(case_data, pop_data, 0.1)
yearly_foi_sd <- 0.6
foi_mean <- logitnorm::twCoefLogitnorm(foi_init(77, median(case_data_indiv$age,na.rm=T), incidence), yearly_foi_sd)[1]
foi_sd <- logitnorm::twCoefLogitnorm(foi_init(77, median(case_data_indiv$age,na.rm=T), incidence), yearly_foi_sd)[2]



# Lower and upper bound of age classes
lr_bound = c(1,seq(5,45,by=5))
ur_bound = c(seq(4,49,by=5))


## Input parameters to model
T_lambda=17 ## How many years of FOI are we trying to estimate? FOI estimates for years before the earliest birth year in the serosurvey will be highly uncertain
n_serotypes=4 ## Assume that n_serotypes of dengue circulate in AS every year. This doesn't really make a difference unless we are trying to estimate the proportion with
## monotypic vs multitypic immunity
covariates=c('sexf') ## What covariates do we want to adjust FOI by?
covariate_names=c('Female sex') ## Names of covariates, for the final results table
svyweights=1 ## Whether to use survey weights
suffix='' ## Suffix to be appended to the file name of the results


## Validation data - number of positive/negative controls, and true positives/false negatives
## From CDC data
N_pos_control = 193
control_tp = 173
N_neg_control = 207
control_fp = 9


maternal_immunity = 3
