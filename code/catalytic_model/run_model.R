## Input data
#
## Serosurvey data
## One line per individual
## Analysis code assumes the following variables:
## yob - the year that the child was born
## birth_month - the month that the child was born (1-12)
## Wgt - survey weight = inverse probability of being selected based on census data
## result_final_no_indt = seropositivity (0/1)
## islandcode = numerical code representing island of residence 
## endres_XXXX where XXXX are years from the earliest birth year in the sample population to the year of the serosurvey
## Values for endres_XXXX are: -1 for years before the child was born; 0 for non-endemic; 1 for endemic outside as; 2 for as
## Note in residence history we do not track the island within as, as this was not provided for all individuals
## covariates: covariates of interest e.g. sex, race. Should be coded as 0/1 for binary, and as dummy variables for categorical (e.g. racewhite, raceblack, raceother)

## Case data and population data

## This case data is derived from a line list of case data with age and year
## The input to the model is a matrix with number of columns equal to the time space of the case data (in this case, 2010-2022) 
## and number of rows equal to the number of age groups considered (in this case, 10 age groups)
## Each entry of case_data is the number of cases in that year and age group

## The population data (pop_data) is a matrix of the same dimension, with the total population of as in each age group and year as each entry
## We assume population is constant over time
## Also require the fraction of the total population that is in each island (i.e. Thomas/John vs Croix)
## We have to assume the age distribution is the same in each island for lack of any data

rm(list=ls())

# setwd("code/catalytic_model/")

source("code/catalytic_model/load_libraries.R")
source("code/catalytic_model/load_case_data.R")
source("code/catalytic_model/load_census_data.R")
source("code/catalytic_model/set_parameters.R")
source("code/catalytic_model/load_serosurvey_data.R")

options(mc.cores = parallel::detectCores()) # checks number of cores without having later to specify the cores argument
rstan_options(auto_write = TRUE) # extended packages to use stan


#### lets try the added frailty model
model_FRAILTYcode = stanc(file="code/catalytic_model/AS_dengue_catalytic_schoolFOI_0513.stan")
compiled_FRAILTYmodel = stan_model(stanc_ret = model_FRAILTYcode)

#run the for each population model

for(pop in names(analysis_set)){

      this_dataset <- analysis_set[[pop]]
      
      as_denguemod = run_FRAILTYregression_model_case_ss(compiled_FRAILTYmodel,
                                                         serosurv_data=this_dataset, #serosurv_data,
                                                         case_data,
                                                         pop_data,
                                                         lr_bound,
                                                         ur_bound,
                                                         svyweights=1,
                                                         c(),
                                                         c(),
                                                         maternal_immunity,#T_lambda,
                                                         n_serotypes,T,glue::glue('_case2016to2022_{pop}_{Sys.Date()}'),
                                                         number_of_chains = 4
                                                         )

      
      lastyear_annualfoi=2023
      serosurv_year=2023
      firstyear_toreport=2007
      
      this_object <- list(
            population = pop,
            dataset= this_dataset,
            fit= as_denguemod,
            p1 = seroprev_fitgraph(this_dataset,as_denguemod),
            p2 = year_fitgraph(as_denguemod),
            p3 = agegroup_fitgraph(as_denguemod),
            foi_table = make_foi_table(as_denguemod,lastyear_annualfoi,serosurv_year,firstyear_toreport)
            
      )
      
      this_object[["combined_plot"]] <- cowplot::plot_grid(this_object$p1,
                                                           cowplot::plot_grid(this_object$p2,this_object$p3),
                                                           ncol=1
                                                           )
      

      write_rds(this_object,file = glue::glue('output/model_fits/{format(Sys.time(), "%Y%m%d")}_{hour(Sys.time())}{minute(Sys.time())}_{pop}.rds'))
            

}



######### combine 2017 and 2018 into one year----
## assume that 2017 never happened
#however, we want to maintain 2023 as the year of the serosurvey
#so increase the year of birth of everyone by 1
#and increase the year of residency for everyone before 2018

#case data
combo_case_data <- case_data[,-2]
combo_case_data[,2] <- case_data[,2] + combo_case_data[,2]

colnames(combo_case_data) <- 2017:2023


#serosurvey data

combo_former_sero <- analysis_set$former %>%
      select(-endres_2017,-res_2017) %>%
      mutate(yob=yob+1)


colnames(combo_former_sero)[2:11] <- paste0("res_",2008:2017)
colnames(combo_former_sero)[33:42] <- paste0("endres_",2008:2017)

T_case <- 7
T_lambda <- 16
as_denguemod = run_FRAILTYregression_model_case_ss(compiled_FRAILTYmodel,
                                                   serosurv_data=combo_former_sero, #modified serosurv_data,
                                                   combo_case_data, #modified case data
                                                   pop_data[,-1], ## modified pop data
                                                   lr_bound,
                                                   ur_bound,
                                                   svyweights=1,
                                                   c(),
                                                   c(),
                                                   maternal_immunity,#T_lambda,
                                                   n_serotypes,T,glue::glue('_case2016to2022_former_combo1718_{Sys.Date()}'),
                                                   number_of_chains = 4,
                                                    combine_17_18=TRUE

)



write_rds(as_denguemod,file = glue::glue('output/model_fits/{format(Sys.time(), "%Y%m%d")}_{hour(Sys.time())}{minute(Sys.time())}_former_combo1718.rds'))



# as_denguemod = run_FRAILTYregression_model_case_ss(compiled_FRAILTYmodel,
#                                                    serosurv_data=this_dataset, #serosurv_data,
#                                                    case_data,
#                                                    pop_data,
#                                                    lr_bound,
#                                                    ur_bound,
#                                                    svyweights=1,
#                                                    c(),
#                                                    c(),
#                                                    maternal_immunity,#T_lambda,
#                                                    n_serotypes,T,glue::glue('_case2016to2022_{pop}_{Sys.Date()}'),
#                                                    number_of_chains = 4
# )



# everyone <- read_rds(file="output/20250507_926_everyone.rds")
# lifelong <- read_rds(file="output/20250507_932_lifelong.rds")
# former <- read_rds(file="output/20250507_938_former.rds")
# latter <- read_rds(file="output/20250507_944_latter.rds")
# 
# 
# everyone$fit %>% rhat() %>% hist()
# lifelong$fit %>% rhat() %>% hist()
# former$fit %>% rhat() %>% hist()
# latter$fit %>% rhat() %>% hist()
# 
# everyone$fit %>% shinystan::launch_shinystan()
# lifelong$fit %>% rhat() %>% hist()
# former$fit %>% rhat() %>% hist()
# latter$fit %>% rhat() %>% hist()
# 
# 
# 
# ## run for 2 to 4 serotypes
# for(circ_sero in 2:4){
#       
#       this_dataset <- analysis_set[["former"]]
#       
#       # as_denguemod = run_FRAILTYregression_model_case_ss(compiled_FRAILTYmodel,
#       #                                                    serosurv_data=this_dataset, #serosurv_data,
#       #                                                    case_data,
#       #                                                    pop_data,
#       #                                                    lr_bound,
#       #                                                    ur_bound,
#       #                                                    svyweights=1,
#       #                                                    c(),
#       #                                                    c(),
#       #                                                    maternal_immunity,T_lambda,
#       #                                                    n_serotypes=circ_sero,T,glue::glue('_case2016to2022_former_{circ_sero}serotypes_{Sys.Date()}'),
#       #                                                    number_of_chains = 4
#       # )
#       
#       as_denguemod = run_FRAILTYregression_model_case_ss(compiled_FRAILTYmodel,
#                                                          serosurv_data=this_dataset, #serosurv_data,
#                                                          case_data,
#                                                          pop_data,
#                                                          lr_bound,
#                                                          ur_bound,
#                                                          svyweights=1,
#                                                          c(),
#                                                          c(),
#                                                          maternal_immunity,#T_lambda,
#                                                          n_serotypes=circ_sero,T,glue::glue('_case2016to2022_former_{circ_sero}serotypes_{Sys.Date()}'),
#                                                          number_of_chains = 4
#       )
#       
#       
#       this_object <- list(
#             population = "former",
#             dataset= this_dataset,
#             fit= as_denguemod,
#             p1 = seroprev_fitgraph(this_dataset,as_denguemod),
#             p2 = year_fitgraph(as_denguemod),
#             p3 = agegroup_fitgraph(as_denguemod),
#             foi_table = make_foi_table(as_denguemod)
#             
#       )
#       
#       this_object[["combined_plot"]] <- cowplot::plot_grid(this_object$p1,
#                                                            cowplot::plot_grid(this_object$p2,this_object$p3),
#                                                            ncol=1
#       )
#       
#       
#       write_rds(this_object,file = glue::glue('output/model_fits/{format(Sys.time(), "%Y%m%d")}_{hour(Sys.time())}{minute(Sys.time())}_former_{circ_sero}serotypes.rds'))
#       
#       
# }





# as_denguemod <- read_rds("as_denguefoi_fixeff_casess_Tlambda_17_mi_3_covs__ns_4_case2016to2022_everyone_2025-02-21.rds")
# 
# as_denguemod
# 
# as_denguemod %>% shinystan::launch_shinystan()
# 
# 
# spread_draws(as_denguemod,frailty[school]) %>% 
#       group_by(school) %>% summarize(median(frailty)) 
# 
# 
# 
# school_labels<- analysis_set$everyone %>% mutate(cluster=as.numeric(as.factor(School)))%>%distinct(School,cluster)
# 
# spread_draws(as_denguemod,frailty[cluster]) %>%
#       left_join(school_labels)%>%
#       ggplot(aes(x=School,group=School,y=frailty))+
#       geom_boxplot()+
#       geom_hline(yintercept = 1,lty=2)+
#       theme(axis.text.x=element_text(angle=45,hjust=1))
#       
#       
# spread_draws(as_denguemod, lambda[i]) %>%
#       mutate(year=2023-57+i) %>%
#       mutate(dataset="2016-2023") %>%
#       filter(year>=2010)%>%
#       ggplot(aes(x=factor(year),y=1-exp(-lambda*4),col=dataset))+
#       geom_boxplot()
# 
# 
# fit_newdata <- readRDS("as_denguefoi_fixeff_casess_Tlambda_17_mi_3_covs__ns_4_case2016to2022.rds")
# fit_newdata_allAS <- readRDS("as_denguefoi_fixeff_casess_Tlambda_17_mi_3_covs__ns_4_case2016to2022_AllAS.rds")
# fit_olddata <- readRDS("as_denguefoi_fixeff_casess_Tlambda_17_mi_3_covs__ns_4_case2010to2022.rds")
# 
# 
# draws_df <-bind_rows(
#       spread_draws(fit_newdata, lambda[i]) %>%
#             mutate(year=2023-57+i) %>%
#             mutate(dataset="2016-2023") ,
#       spread_draws(fit_olddata, lambda[i]) %>%
#             mutate(year=2023-63+i)%>%
#             mutate(dataset="2010-2023"),
#       spread_draws(fit_newdata_allAS, lambda[i]) %>%
#             mutate(year=2023-57+i) %>%
#             mutate(dataset="2016-2023 (All AS)") 
#       )
# 
# 
# draws_df %>%
#       filter(year>=2010)%>%
#       ggplot(aes(x=factor(year),y=lambda*4,col=dataset))+
#       geom_boxplot()
# 
# spread_draws(fit_olddata, lambda_temp[i]) %>%
#       mutate(year=1961 + i-1)%>%
#       mutate(dataset="2010-2023") %>%
#       group_by(year,i)%>%
#       summarize(mean(4*lambda_temp))%>%
#       arrange(-year) %>%
#       View()
      

# model_code = stanc(file="AS_dengue_catalytic.stan")
# compiled_model = stan_model(stanc_ret = model_code)
# 
# #rerun the model
# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                             serosurv_data,
#                                             case_data,
#                                             pop_data,
#                                             lr_bound,
#                                             ur_bound,
#                                             svyweights=1,
#                                             c(),
#                                             c(),
#                                             maternal_immunity,T_lambda,
#                                             n_serotypes,T,'_case2010to2022')
# 
# 
# # write_rds(as_denguemod,"as_denguemod.rds")
# 
# as_denguemod %>%
#       shinystan::launch_shinystan()
# 
# 
# #use the saved model
# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                             serosurv_data,
#                                             case_data,
#                                             pop_data,
#                                             # totalpop=sum(popdata$npop),
#                                             lr_bound,
#                                             ur_bound,
#                                             svyweights=1,
#                                             c(),
#                                             c(),
#                                             maternal_immunity,T_lambda,
#                                             n_serotypes,F,'_case2010to2022')





## OLD




#############################
## Analysis of dengue serosurvey in as
## Matt Hitchings 09/08/22
# require(dplyr)
# require(rstan)
# require(readxl)
# require(excel.link)
# require(tidyverse)
# require(bayesplot)
# require(binom)
# require(data.table)
# require(grid)
# require(gridExtra)
# require(ggpubr)

## Input data

## Serosurvey data
## One line per individual
## Analysis code assumes the following variables:
## yob - the year that the child was born
## birth_month - the month that the child was born (1-12)
## Wgt - survey weight = inverse probability of being selected based on census data
## result_final_no_indt = seropositivity (0/1)
## islandcode = numerical code representing island of residence 
## endres_XXXX where XXXX are years from the earliest birth year in the sample population to the year of the serosurvey
## Values for endres_XXXX are: -1 for years before the child was born; 0 for non-endemic; 1 for endemic outside as; 2 for as
## Note in residence history we do not track the island within as, as this was not provided for all individuals
## covariates: covariates of interest e.g. sex, race. Should be coded as 0/1 for binary, and as dummy variables for categorical (e.g. racewhite, raceblack, raceother)

## Case data and population data

## This case data is derived from a line list of case data with age and year
## The input to the model is a matrix with number of columns equal to the time space of the case data (in this case, 2010-2022) 
## and number of rows equal to the number of age groups considered (in this case, 10 age groups)
## Each entry of case_data is the number of cases in that year and age group

## The population data (pop_data) is a matrix of the same dimension, with the total population of as in each age group and year as each entry
## We assume population is constant over time
## Also require the fraction of the total population that is in each island (i.e. Thomas/John vs Croix)
## We have to assume the age distribution is the same in each island for lack of any data

# # Case data
# case_data_indiv = xl.read.file("Data/ArboNET Case Data.xlsx")
# case_data_indiv$age = as.numeric(case_data_indiv$age)
# 
# cases_by_year = case_data_indiv %>% group_by(year) %>% summarise(ncase = n())
# # cases_by_year = rbind(cases_by_year,data.frame('year'=c(2010:2015,2019:2021,2023),'ncase'=0)) %>% arrange(year)
# #
# cases_by_year = rbind(cases_by_year,data.frame('year'=c(2019:2021,2023),'ncase'=0)) %>% arrange(year)
# 
# 
# 
# allagegroups = paste0(seq(0,80,by=5)," to ",seq(4,84,by=5))
# agegroups = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
# 
# caseage_by_year = case_data_indiv %>% filter(age<50 & age>0) %>%
#   mutate(AgeGroup = cut(age,breaks=seq(-0.01,50,by=5),labels=agegroups)) %>% 
#   group_by(year,AgeGroup) %>% summarise(ncase=n())
# 
# AG=10 # Number of age groups considered
# max_A=49 # Max age considered in the case (note, does not have to be max age. The oldest cases may provide little info about transmission, as we
# ## are primarily interested in FOI over the period that serosurvey participants were alive)
# # T_case=14 # Time span of case data
# #
# T_case=8 # Time span of case data



# # reporting parameter priors
# reporting_mean <- -2.2
# reporting_sd <- 2

# # Function to estimate FOI based on mean age of infection
# foi_init <- function(L, A, inc) {
#   r0 <- 1 + L / A
#   s <- 1 / r0
#   inc / s
# }
# 
# # Function to estimate incidence rate from case data
# calc_incidence <- function(cases, pop, reporting_rate) {
#   if (is.null(dim(pop))) {
#     N <- sum(pop) # pop size
#   } else {
#     N <- sum(pop[, 1]) # pop size
#   }
#   total_cases <- sum(cases)
#   nb_year <- ncol(cases) # years available of incidence data
#   inc <- total_cases / nb_year / N / reporting_rate
#   return(inc)
# }

# # Construct case data to input to model
# case_data = matrix(nrow=AG,ncol=T_case)
# agegroups = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
# rownames(case_data) = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
# # colnames(case_data) = 2010:2023
# #
# colnames(case_data) = 2016:2023
# 
# 
# for (r in 1:nrow(case_data)) {
#   for (c in 1:ncol(case_data)) {
#     case_data[r,c] = ifelse(length(caseage_by_year$ncase[caseage_by_year$year==colnames(case_data)[c] & caseage_by_year$AgeGroup==rownames(case_data)[r]])>0,
#                             caseage_by_year$ncase[caseage_by_year$year==colnames(case_data)[c] & caseage_by_year$AgeGroup==rownames(case_data)[r]],0)
#   }
# }

# ## Population data comes from US census data
# popdata = read_excel('Data/2020 census.xlsx')
# popdata = popdata %>% 
#   mutate(AgeGroup = gsub(' years','',AgeGroup)) %>% 
#   filter(AgeGroup %in% agegroups) %>% mutate(npop=Number)
# 
# pop_data = matrix(nrow=AG,ncol=T_case)
# rownames(pop_data) = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
# # colnames(case_data) = 2010:2023
# #
# colnames(pop_data) = 2016:2023
# for (r in 1:nrow(pop_data)) {
#   for (c in 1:ncol(pop_data)) {
#     pop_data[r,c] = popdata$npop[popdata$AgeGroup==rownames(pop_data)[r]]
#   }
# }



# ## Code to get prior parameters for the yearly FOI, based on observed incidence
# incidence <- calc_incidence(case_data, pop_data, 0.1)
# yearly_foi_sd <- 0.6
# foi_mean <- logitnorm::twCoefLogitnorm(foi_init(77, median(case_data_indiv$age,na.rm=T), incidence), yearly_foi_sd)[1]
# foi_sd <- logitnorm::twCoefLogitnorm(foi_init(77, median(case_data_indiv$age,na.rm=T), incidence), yearly_foi_sd)[2]
# 
# 
# 
# # Lower and upper bound of age classes
# lr_bound = c(1,seq(5,45,by=5))
# ur_bound = c(seq(4,49,by=5))



# ## Serosurvey data
# serosurv_data = read.csv('Data/merged_latter1.csv') %>% dplyr::select(-X) %>% rename(RecordID=Record.ID,
#                                                                                        Age=Calculated.age.based.on.DOB.,
#                                                                                        Village=Village.of.Residence) %>%
#       filter(Final.Test.Result!="I")
# colnames(serosurv_data) = gsub('X','res_',colnames(serosurv_data))
# 
# survweights = read.csv('Data/AS Data_Probability.csv')
# 
# serosurv_data = serosurv_data %>% left_join(survweights,by='RecordID')
# serosurv_data = serosurv_data %>% mutate(Seropositive = as.numeric(Final.Test.Result=="P"),
#                                          sexf = as.numeric(Sex=="Female"))
# 
# ## Define endemicity of residence by year
# allcountries = unique(unlist(c(serosurv_data[,paste0('res_',2007:2022)])))
# endemic_list = c("Fiji","Philippines","Western Samoa","Tuvalu")
# 
# for (y in 2007:2023) {
#   
#   serosurv_data[,paste0("endres_",y)] = -1*as.numeric(is.na(serosurv_data[,paste0("res_",y)])) + 
#     as.numeric(serosurv_data[,paste0("res_",y)] %in% endemic_list) + 
#     2*as.numeric(serosurv_data[,paste0("res_",y)] %in% c("American Samoa"))
#   
# }


# ## Input parameters to model
# T_lambda=17 ## How many years of FOI are we trying to estimate? FOI estimates for years before the earliest birth year in the serosurvey will be highly uncertain
# n_serotypes=4 ## Assume that n_serotypes of dengue circulate in AS every year. This doesn't really make a difference unless we are trying to estimate the proportion with
# ## monotypic vs multitypic immunity
# covariates=c('sexf') ## What covariates do we want to adjust FOI by?
# covariate_names=c('Female sex') ## Names of covariates, for the final results table
# svyweights=1 ## Whether to use survey weights
# suffix='' ## Suffix to be appended to the file name of the results
# 
# 
# ## Validation data - number of positive/negative controls, and true positives/false negatives
# ## From CDC data
# N_pos_control = 193
# control_tp = 173
# N_neg_control = 207
# control_fp = 9

## How many months infants are assumed to be immune to dengue after they are born (to adjust FOI in their first year of life)
## I found this had little impact on results
# maternal_immunity = 3


# # Observed seroprevalence and binomial CIs by birth year and island
# bb_data = serosurv_data %>% group_by(yob) %>% summarise(n=n(),x=sum(Seropositive))
# bin_cis = t(sapply(1:nrow(bb_data),
#                    function(x) {
#                      d=binom.confint(bb_data$x[x],bb_data$n[x])
#                      return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
#                    }))
# colnames(bin_cis)=c("SP","lwr","upr")

##### Function to run model
# expit = function(x) {1/(1+exp(-x))}
# run_regression_model_case_ss = function(compiled_model,
#                                         serosurv_data,
#                                         case_data,
#                                         pop_data,
#                                         lr_bound,
#                                         ur_bound,
#                                         svyweights,
#                                         covariates,covariate_names,
#                                         maternal_immunity,T_lambda,
#                                         n_serotypes,runmod,suffix,varybyisland=T) {
# 
#   # Define weights based on whether we are using survey weights or not
#   serosurv_data$nonw = 1
#   if (svyweights==1) {
#     weights = 1/serosurv_data[,'prob']
#   } else {
#     weights = serosurv_data[,'nonw']
#   }
# 
# 
#   # Pick out covariates from serosurvey data based on names
#   if (length(covariates)>0) {
# 
#         print(length(covariates))
#     reg_vars=serosurv_data[,covariates,drop=FALSE]
#   } else {
#     reg_vars=array(0,dim=c(nrow(serosurv_data),0))
#   }
#   r_v = length(covariates)
# 
#   # Unique file suffix based on input parameters
#   filesuffix = paste0('_fixeff_casess_Tlambda_',T_lambda,'_mi_',maternal_immunity,
#                       '_covs_',paste0(covariates,collapse='.'),
#                       '_ns_',n_serotypes,
#                       suffix)
# 
#   print(filesuffix)
# 
#   # Input data for rstan
#   data_prov <-
#     list(
#       N = nrow(serosurv_data), # the number of individuals in the serosurvey
#       A_sp = 2023 - min(serosurv_data$yob) + 1, #number of age classes in the serosurvey
#       Y = serosurv_data$Seropositive, #serostatus of each individual
#       w = weights,
#       endres_history = as.matrix(serosurv_data[,paste0("endres_",2007:2023)]),
#       r_v = r_v,
#       reg_vars = reg_vars,
#       birthyear = serosurv_data$yob-2007+1, # yob goes in as number of years since the earliest birth year
#       birthmonth = month(serosurv_data$Date.of.birth),
#       T_ss = T_lambda, #number of years for which we are estimating FOI
#       N_pos_control = N_pos_control,
#       control_tp = control_tp,
#       N_neg_control = N_neg_control,
#       control_fp = control_fp,
#       maternal_immunity = maternal_immunity,
#       n_serotypes = n_serotypes,
#       AG=AG,
#       max_A=max_A,
#       T_case=T_case,
#       secondary_cases=case_data,
#       pop=pop_data,
#       lr_bound=lr_bound,
#       ur_bound=ur_bound,
#       foi_mean=foi_mean,
#       foi_sd=foi_sd,
#       reporting_mean=reporting_mean,
#       reporting_sd=reporting_sd,
#       pos2010_18to25=179,
#       total2010_18to25 = 201,
#       pos2010_26to40 = 216,
#       total2010_26to40 = 217
# 
# 
#          )
# 
#   ### Without starting values
#   if (runmod) {
#     system.time(
#       as_denguefoi_fe <-
#         sampling(
#           compiled_model,
#           data = data_prov,
#           chains = 4,
#           iter = 1000,
#           warmup = 500,
#           thin = 1,
#           control = list(adapt_delta = 0.80)#,
#           # init = list(
#           #       list(exp_sp_18to25=1,exp_sp_26to40=1),
#           #       list(exp_sp_18to25=1,exp_sp_26to40=1),
#           #       list(exp_sp_18to25=1,exp_sp_26to40=1),
#           #       list(exp_sp_18to25=1,exp_sp_26to40=1)
#           # )
# 
#           )
#     )
#     as_denguefoi_fe@stanmodel@dso <- new("cxxdso")
#     saveRDS(as_denguefoi_fe, file = paste0('as_denguefoi',filesuffix,'.rds'))
# 
#   } else {
#     as_denguefoi_fe = readRDS(paste0('as_denguefoi',filesuffix,'.rds'))
#   }
# 
#   ## The rest of the code is processing, making plots, etc.
#   ## Without the variation by island, this code will have to be changed. I can work on this after
#   ## you are more familiar with the code
# 
#   exp_sp = rstan::extract(as_denguefoi_fe,pars="exp_sp", inc_warmup=FALSE, permute=T)$exp_sp
#   exp_sp_mean <- apply(exp_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
# 
#   exp_obs_sp = rstan::extract(as_denguefoi_fe,pars="exp_obs_sp", inc_warmup=FALSE, permute=T)$exp_obs_sp
#   exp_obs_sp_mean <- apply(exp_obs_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
# 
#   exp_obs_sp_byage = sapply(1:2000,function(i) {
#     serosurv_data$exp_obs_sptemp=exp_obs_sp[i,]
#     exp_obs_sp_byage = serosurv_data %>% group_by(yob) %>% summarise(esp=mean(exp_obs_sptemp))
#     return(exp_obs_sp_byage$esp)
#   })
#   exp_obs_sp_byage_summ = apply(exp_obs_sp_byage,1,function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   lambda = rstan::extract(as_denguefoi_fe,pars="lambda_temp",inc_warmup=F,permute=T)$lambda_temp
#   lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
#   colnames(lambda)=1961 + (0:(ncol(lambda)-1))
# 
#   cum_lambda_sp = exp(rstan::extract(as_denguefoi_fe,pars="log_cum_lambda",inc_warmup=F,permute=F))
#   cum_lambda_sp = apply(cum_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
#   end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
#   nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   bb_data = serosurv_data %>% group_by(yob) %>% summarise(n=n(),x=sum(Seropositive))
#   bin_cis = as.data.frame(t(sapply(1:nrow(bb_data),
#                      function(x) {
#                        d=binom.confint(bb_data$x[x],bb_data$n[x])
#                        return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
#                      })))
#   colnames(bin_cis)=c("SP","lwr","upr")
#   bin_cis$yob = bb_data$yob
# 
#   d_forplot = data.frame('Age'=rep(bin_cis$yob,2),
#                          'Model'=rep(c('Observed','Model estimated'),each=nrow(bin_cis)),
#                          'SP'=c(bin_cis[,1],
#                                 exp_obs_sp_byage_summ[2,]),
#                          'lwr'=c(bin_cis[,2],exp_obs_sp_byage_summ[1,]),
#                          'upr'=c(bin_cis[,3],exp_obs_sp_byage_summ[3,]))
# 
#   theme_set(theme_classic())
#   psero=ggplot() +
#     #geom_errorbar(data = d_forplot,aes(x=2022-Age+1,ymin=lwr,ymax=upr,col=factor(Model,levels=c('Observed','Model estimated'))),width=0) +
#     geom_line(data = d_forplot,aes(x=2023-Age+1,y=SP,col=factor(Model,levels=c('Observed','Model estimated'))))+
#     geom_ribbon(data = d_forplot,aes(x=2023-Age+1,ymin=lwr,ymax=upr,fill=factor(Model,levels=c('Observed','Model estimated'))),alpha=0.3)+
#     geom_point(data = d_forplot, aes(x=2023-Age+1,y=SP,shape=factor(Model,levels=c('Observed','Model estimated'))))+
#     scale_color_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
#     scale_fill_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
#     scale_shape_manual(name=NULL,values=c(16,NA),labels=c('Observed','Model estimated'))+
#     ylab('Seroprevalence')+ylim(c(0,1))+
#     scale_x_continuous(name='Birth year',breaks=9:17,labels=2015:2007)+
#     theme(axis.text=element_text(size=6),
#           axis.title = element_text(size=6),
#           legend.title = element_text(size=6),
#           legend.text = element_text(size=6),
#           legend.position = c(0.1,0.89),
#           strip.text = element_text(size=6))
# 
#   ggsave(paste0('./ModelFitBirthYear',filesuffix,'.png'),
#          psero,
#          device='png',height=6,width=6)
# 
#   reprate = rstan::extract(as_denguefoi_fe,'reporting_rate',inc_warmup=F,permute=T)$reporting_rate
#   reprate = apply(reprate,c(2),function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   susc = rstan::extract(as_denguefoi_fe,'susc',inc_warmup=F,permute=T)$susc
#   susc = apply(susc,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   mono = rstan::extract(as_denguefoi_fe,'mono',inc_warmup=F,permute=T)$mono
#   mono = apply(mono,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
# 
#   repcases = rstan::extract(as_denguefoi_fe,'exp_reported_cases')$exp_reported_cases
#   repcases = apply(repcases,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
#   repcases_forplot = as.data.frame(apply(repcases,1,cbind))
#   colnames(repcases_forplot)=c("LCI","Median","UCI")
#   repcases_forplot$AgeGroup=agegroups
#   # repcases_forplot$Year=rep(2010:2023,each=10)
#   repcases_forplot$Year=rep(2016:2023,each=10)
# 
# 
# 
#   repcasestotal_byyear = apply(repcases,c(1,3),sum)
#   repcasestotal_byyear = as.data.frame(t(apply(repcasestotal_byyear,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
#   colnames(repcasestotal_byyear)=c("LCI","Median","UCI")
#   # repcasestotal_byyear$Year=2010:2023
#   repcasestotal_byyear$Year=2016:2023
# 
#   repcasestotal_byyear$Model='Estimated'
# 
#   obscasestotal = data.frame('Median'=colSums(case_data),
#                              'LCI'=NA,
#                              'UCI'=NA,
#                              # 'Year'=2010:2023,
#                              'Year'=2016:2023,
#                              'Model'="Observed")
#   repcasestotal_byyear = rbind(repcasestotal_byyear,obscasestotal)
# 
# 
#   p1=ggplot() + geom_col(data=repcasestotal_byyear %>% filter(Model=="Observed"),aes(x=Year,y=Median),col='grey',alpha=.3)+
#     geom_point(data=repcasestotal_byyear %>% filter(Model=="Estimated"),aes(x = Year, y = Median),col = "steelblue") +
#     geom_errorbar(data=repcasestotal_byyear %>% filter(Model=="Estimated"),
#                 aes(x=Year,ymin=LCI,ymax=UCI),col = "steelblue")+
#     # scale_x_continuous(name="Year",breaks=seq(2010,2023,by=2))+ylab("Reported cases")+
#         scale_x_continuous(name="Year",breaks=seq(2016,2023,by=2))+ylab("Reported cases")+
#     theme(axis.text=element_text(size=6),
#           axis.title=element_text(size=6))
# 
#   repcasestotal_byage = apply(repcases,c(1,2),sum)
#   repcasestotal_byage = as.data.frame(t(apply(repcasestotal_byage,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
#   colnames(repcasestotal_byage)=c("LCI","Median","UCI")
#   repcasestotal_byage$AgeGroup=agegroups
#   repcasestotal_byage$Model='Estimated'
# 
#   obscasestotal = data.frame('Median'=rowSums(case_data),
#                              'LCI'=NA,
#                              'UCI'=NA,
#                              'AgeGroup'=agegroups,
#                              'Model'="Observed")
#   repcasestotal_byage = rbind(repcasestotal_byage,obscasestotal)
# 
#   repcasestotal_byage$AgeGroup = factor(repcasestotal_byage$AgeGroup,
#                                         levels=agegroups)
# 
#   agegroups_plot = c("1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49")
#   p2=ggplot() + geom_col(data=repcasestotal_byage %>% filter(Model=="Observed"),aes(x=AgeGroup,y=Median),col='grey',alpha=.3)+
#     geom_point(data=repcasestotal_byage %>% filter(Model=="Estimated"),aes(x = AgeGroup, y = Median),col = "steelblue") +
#     geom_errorbar(data=repcasestotal_byage %>% filter(Model=="Estimated"),
#                   aes(x=AgeGroup,ymin=LCI,ymax=UCI),col = "steelblue")+
#     scale_x_discrete(name="Age Group",labels=agegroups_plot,guide = guide_axis(n.dodge = 2))+
#     ylab("Reported cases")+
#     theme(axis.text=element_text(size=6),
#           axis.title=element_text(size=6))
#   ggsave(paste0('./ModelFitSeroprev_Cases',filesuffix,'.png'),
#          grid.arrange(psero,grid.arrange(p1,p2,ncol=2),nrow=2,heights=c(4,2)),
#          device='png',height=6,width=6,units='in')
# 
#   # ggplot(repcases_forplot, aes(x = AgeGroup, y = Median), col = "steelblue", alpha = .3) +
#   #   geom_col() +
#   #   facet_wrap(vars(Year))+
#   #   theme_classic()
# 
#   fitarray=as.array(as_denguefoi_fe)
# 
#   fitarray[,,grepl('lambda_temp',dimnames(fitarray)[[3]])] =
#     n_serotypes*(fitarray[,,grepl('lambda_temp',dimnames(fitarray)[[3]])])
#   fitarray[,,grepl('beta_covs',dimnames(fitarray)[[3]])] = exp(fitarray[,,grepl('beta_covs',dimnames(fitarray)[[3]])])
#   fitarray[,,'log_cum_lambda'] = n_serotypes*exp(fitarray[,,'log_cum_lambda'])
#   fitarray[,,'logit_lambda_end'] = n_serotypes*expit(fitarray[,,'logit_lambda_end'])
#   fitarray[,,'logit_lambda_nonend'] = n_serotypes*expit(fitarray[,,'logit_lambda_nonend'])
# 
# 
#   print(r_v)
# 
#   # if (r_v>0) {
#   #   if (data_prov$T_ss<data_prov$A_sp) {
#   #     fitarray_forplot = fitarray[,,c(paste0('lambda_temp[',48:59,']'),'log_cum_lambda',
#   #                                     'hr[1]',
#   #                                     'logit_lambda_end','logit_lambda_nonend',
#   #                                     'spec','sens',paste0('beta_covs[',1:r_v,']')
#   #     )]
#   #     dimnames(fitarray_forplot)[[3]]=c(2008:2019,'Cumulative FOI 2020-2022',
#   #                                       'HR (St John/St Thomas vs St Croix)',
#   #                                       'Endemic','Non-endemic','Specificity','Sensitivity'
#   #                                       ,paste0('HR: ',covariate_names)
#   #     )
#   #   } else {
#   #     fitarray_forplot = fitarray[,,c(paste0('lambda_temp[',48:62,']'),
#   #                                     'hr[1]',
#   #                                     'logit_lambda_end','logit_lambda_nonend',
#   #                                     'spec','sens',paste0('beta_covs[',1:r_v,']')
#   #     )]
#   #     dimnames(fitarray_forplot)[[3]]=c(2008:2022,
#   #                                       'HR (St John/St Thomas vs St Croix)',
#   #                                       'Endemic','Non-endemic','Specificity','Sensitivity'
#   #                                       ,paste0('HR: ',covariate_names)
#   #     )
#   #   }
#   #
#   # } else {
#   #   if (data_prov$T_ss<data_prov$A_sp) {
#   #     fitarray_forplot = fitarray[,,c(paste0('lambda_temp[',48:59,']'),'log_cum_lambda',
#   #                                     'hr[1]',
#   #                                     'logit_lambda_end','logit_lambda_nonend',
#   #                                     'spec','sens'
#   #     )]
#   #     dimnames(fitarray_forplot)[[3]]=c(2008:2019,'Cumulative FOI 2020-2022',
#   #                                       'HR (St John/St Thomas vs St Croix)',
#   #                                       'Endemic','Non-endemic','Specificity','Sensitivity')
#   #   } else {
#   #     fitarray_forplot = fitarray[,,c(paste0('lambda_temp[',48:62,']'),
#   #                                     'hr[1]',
#   #                                     'logit_lambda_end','logit_lambda_nonend',
#   #                                     'spec','sens'
#   #     )]
#   #     dimnames(fitarray_forplot)[[3]]=c(2008:2022,
#   #                                       'HR (St John/St Thomas vs St Croix)',
#   #                                       'Endemic','Non-endemic','Specificity','Sensitivity')
#   #   }
#   # }
#   #
#   # theme_set(theme_classic())
#   # ggsave(paste0('./plots/MCMCDens',filesuffix,'.png'),
#   #        mcmc_dens(fitarray_forplot)+facet_text(size=6),
#   #        height=7,width=7,units='in',device='png')
# 
#   foi_table = data.frame('Year'=c(2007:(2007+data_prov$T_ss-1),
#                                   rep('Cumulative 2020-2022',as.numeric(data_prov$T_ss<data_prov$A_sp)),
#                                   'Endemic','Non-endemic'))
#   foi_table$FOI = c(paste0(round(n_serotypes*lambda[2,46:(46+data_prov$T_ss-1)],2),' (',
#                            round(n_serotypes*lambda[1,46:(46+data_prov$T_ss-1)],2),',',
#                            round(n_serotypes*lambda[3,46:(46+data_prov$T_ss-1)],2),')'),
#                     rep(paste0(round(n_serotypes*cum_lambda_sp[2],2),' (',
#                            round(n_serotypes*cum_lambda_sp[1],2),',',
#                            round(n_serotypes*cum_lambda_sp[3],2),')'),as.numeric(data_prov$T_ss<data_prov$A_sp)),
#                     paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
#                            round(n_serotypes*end_lambda_sp[1],2),',',
#                            round(n_serotypes*end_lambda_sp[3],2),')'),
#                     paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
#                            round(n_serotypes*nonend_lambda_sp[1],2),',',
#                            round(n_serotypes*nonend_lambda_sp[3],2),')')
#   )
# 
# 
# 
#   write.csv(foi_table,paste0('./FOIByYear',filesuffix,'.csv'),row.names=F)
# 
# 
#   # Regression coefficients
#   reg_coeffs = rstan::extract(as_denguefoi_fe,'beta_covs')$beta_covs
# 
#   if (r_v>0) {
#     rv_table = data.frame('Variable'=c(covariate_names),
#                           'HR_CrI' = apply(reg_coeffs,2,function(x)
#                             paste0(sprintf('%.2f',median(exp(x)))," (",
#                                    sprintf('%.2f',quantile(exp(x),0.025)),",",
#                                    sprintf('%.2f',quantile(exp(x),0.975)),")")
#                           )
#     )
#     write.csv(rv_table,paste0('./DengueasHRs',filesuffix,'.csv'),row.names=F)
#   } else {
# 
#   }
# 
#   return(as_denguefoi_fe)
# 
# }

# model_code = stanc(file="AS_dengue_catalytic.stan")
# compiled_model = stan_model(stanc_ret = model_code)
# 
# #rerun the model
# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                             serosurv_data,
#                                             case_data,
#                                             pop_data,
#                                             # totalpop=sum(popdata$npop),
#                                             lr_bound,
#                                             ur_bound,
#                                             svyweights=1,
#                                             c(),
#                                             c(),
#                                             maternal_immunity,T_lambda,
#                                             n_serotypes,T,'_case2016to2022')
# 
# 
# # write_rds(as_denguemod,"as_denguemod.rds")
# 
# as_denguemod %>%
#       shinystan::launch_shinystan()
# 
# 
# #use the saved model
# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                               serosurv_data,
#                                               case_data,
#                                               pop_data,
#                                               # totalpop=sum(popdata$npop),
#                                               lr_bound,
#                                               ur_bound,
#                                               svyweights=1,
#                                               c(),
#                                               c(),
#                                               maternal_immunity,T_lambda,
#                                               n_serotypes,F,'_case2016to2022')
# library(tidybayes)
# 
# as_denguemod %>% 
#       spread_draws(reporting_rate[t]) %>% 
#       group_by(t)%>%
#       median_qi() %>%
#       mutate(year=2015+t)%>%
#       ggplot(aes(x=year,y=reporting_rate))+
#       geom_ribbon(aes(x=year,ymin = .lower,ymax=.upper),
#                   alpha=0.1
#       )+
#       geom_point()+
#       scale_y_continuous(trans="log10")
# 
# 
# as_denguemod %>% 
#       spread_draws(lambda[t]) %>% 
# group_by(t)%>%
#       mutate(year=2023-57+t)%>%
#       filter(year>=2010)%>%
#       median_qi() %>%
# mutate(lambda=lambda*4)#%>%
# #       pull(lambda)%>%
# #       plot()
# 
# 
# as_denguemod %>% 
#       spread_draws(lambda[t]) %>% 
#       mutate(year=2023-57+t)%>%
#       filter(year>=2011) %>%
#       mutate(period=case_when(
#             year<2016 ~ 1,
#             year %in% c(2016,2017,2018) ~ 2,
#             year>2018 ~ 3
#       ))%>%
#       # mutate(period_yrs=case_when(
#       #       year<2017 ~ 6,
#       #       year %in% c(2017,2018) ~ 3,
#       #       year>2018 ~ 5
#       # ))%>%
#       group_by(.draw,period) %>%
#       summarize(prob=1-prod(1-4*lambda))  %>%
#       group_by(period)%>%
#       summarize(
#                 median=median(prob),
#                 lb=quantile(prob,0.025),
#                 ub=quantile(prob,0.975)
#       ) %>%
#       mutate(period=factor(period,labels=c("Before (2011-2015)",
#                                            "During (2016-2018)",
#                                            "After (2019-2023)"
#       )))%>%
#       ggplot(aes(x=period))+
#       geom_point(aes(y=median))+
#       geom_linerange(aes(ymin=lb,ymax=ub))+
#       scale_y_continuous("Average Annual FOI",limits=c(0,0.5))+
#       xlab("Period relative to outbreak")


### OLD OLD 

# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                               serosurv_data,
#                                               case_data,
#                                               pop_data,
#                                               # totalpop,
#                                               lr_bound,
#                                               ur_bound,
#                                               svyweights,
#                                               c('sexf'),
#                                               c('Female sex'),
#                                               maternal_immunity,T_lambda,
#                                               n_serotypes,T,'_case2010to2022')
# 
# # Sensitivity analyses
# as_denguemod = run_regression_model_case_ss(compiled_model,
#                                               serosurv_data %>% 
#                                                 filter(residence_as=="all_life") %>% 
#                                                 mutate(racewhite = as.numeric(race_cat_red=='white'),
#                                                        raceother = as.numeric(race_cat_red=='other')),
#                                               case_data,
#                                               pop_data,
#                                               totalpop,
#                                               lr_bound,
#                                               ur_bound,
#                                               svyweights,
#                                               c('sexf'),
#                                               c('Female sex'),
#                                               maternal_immunity,T_lambda,
#                                               n_serotypes,T,'_residents_cases2010to2022')
# 
# 
# 
# ### FOI by island
# filesuffix = '_fixeff_casess_Tlambda_15_mi_3_covs_sexf.racewhite.raceother.ethnicity_hi.sch_private_ns_4_case2010to2022'
# as_denguemod = readRDS('as_denguefoi_fixeff_casess_Tlambda_15_mi_3_covs_sexf.racewhite.raceother.ethnicity_hi.sch_private_ns_4_case2010to2022.rds')
# 
# lambda = rstan::extract(as_denguemod,pars="lambda_temp",inc_warmup=F,permute=T)$lambda_temp
# lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
# colnames(lambda)=1961 + (0:(ncol(lambda)-1))
# lambda = lambda[,50:62]
# 
# 
# lambda_byis = rstan::extract(as_denguemod,'lambda_rep',inc_warmup=F,permute=T)$lambda_rep
# lambda_byis = apply(lambda_byis,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
# dimnames(lambda_byis) = list(c('2.5%','50%','97.5%'),2008:2022,c('St Croix','St John/St Thomas'))
# lambda_byis = lambda_byis[,3:15,]
# 
# foi_table_byis = data.frame('Year'=c(2010:2022))
# foi_table_byis$FOI_StCroix = c(paste0(round(n_serotypes*lambda_byis[2,,1],2),' (',
#                          round(n_serotypes*lambda_byis[1,,1],2),',',
#                          round(n_serotypes*lambda_byis[3,,1],2),')')
# )
# foi_table_byis$FOI_StJT = c(paste0(round(n_serotypes*lambda_byis[2,,2],2),' (',
#                                       round(n_serotypes*lambda_byis[1,,2],2),',',
#                                       round(n_serotypes*lambda_byis[3,,2],2),')')
# )
# 
# write.csv(foi_table_byis,paste0('./ParamTables/FOIByYearByIsland',filesuffix,'.csv'),row.names=F)
# 
# 
# 
# 
# # FOI plot
# d_forfoiplot = bind_rows(data.frame(t(n_serotypes*lambda)),data.frame(t(n_serotypes*lambda_byis[,,1])),data.frame(t(n_serotypes*lambda_byis[,,2])))
# colnames(d_forfoiplot) = c('lci','med','uci')
# rownames(d_forfoiplot) = NULL
# d_forfoiplot = d_forfoiplot %>% 
#   mutate(year = rep(2010:2022,3),
#          island = rep(c('as','St. Croix','St. Thomas/St. John'),each=13),
#          repcases = c(c(0,0,cases_by_year$ncase,1,2,0),
#                       colSums(case_data_byis[,,1]),
#                       colSums(case_data_byis[,,2]))
#   )
# theme_set(theme_bw())
# ggsave('./WriteUp/ReportedCasesAndFOIByIsland_1115.png',
#        d_forfoiplot %>% 
#   ggplot() + geom_line(aes(x=year,y=med,color=factor(island)),show.legend = F) + 
#   geom_ribbon(aes(x=year,ymin=lci,ymax=uci,fill=factor(island),color=NULL),alpha=0.5,show.legend = F)+
#   geom_line(aes(x=year,y=repcases/200))+
#   facet_wrap(vars(factor(island,levels=c('as','St. Croix','St. Thomas/St. John'))),ncol=3)+
#   scale_x_continuous(name='Year',breaks=seq(2010,2022,by=4),labels=c('\'10','\'14','\'18','\'22'))+
#   scale_y_continuous(name='Force of infection',limits=c(0,1),sec.axis = sec_axis( trans=~.*200, name="Reported cases")),
#   height=3,width=6,units='in',device='png')
# 
# invroot_phi = 5
# max_A = data_prov$max_A
# T_ss = data_prov$T_ss
# AG = data_prov$AG
# A_sp = data_prov$A_sp
# 
# reporting_rate_logit = rep(NA,T_case)
# lambdaRE_logit = rep(NA,max_A+T_case)
# 
# alpha = 0.3
# 
# reporting_rate0 = rnorm(1,reporting_mean, sqrt(reporting_sd^2 / 2))
# for (t in 1:T_case) {
#   reporting_rate_logit[t] = rnorm(1,reporting_rate0, sqrt(reporting_sd^2 / 2))
# }
# 
# lambda0_logit = rnorm(1,foi_mean, sqrt(foi_sd ^2 / 2))
# for (i in 1:(max_A+T_case)) {
#   lambdaRE_logit[i] = rnorm(1,lambda0_logit, sqrt(foi_sd ^2 / 2))
# }
# 
# 
# phi = 1/(invroot_phi)^2
# 
# reporting_rate = expit(reporting_rate_logit)
# 
# 
# lambda = expit(lambdaRE_logit)
# 
# 
# 
# # the rest of this code is focused around constructing individual-level cumulative FOI based on each 
# # individual's age and residence history
# lambda_i = rep(NA,data_prov$N)# Lambda for each individual included in the serosurvey
# exp_sp = rep(NA,data_prov$N)# Probability of seropositivity for each individual included in the serosurvey
# 
# lambda_rep = rep(NA,T_ss) # Population-representative FOI for AS for time period of serosurvey
# 
# lambda_temp = rep(NA,T_case+max_A) # Temporary to store lambda/lambda_rep
# 
# cum_lambda = matrix(NA,nrow=max_A, ncol=T_case) # Cumulative lambda experienced by residents of AS
# mono = matrix(NA,nrow=max_A, ncol=T_case) # Proportion monotypic
# susc = matrix(NA,nrow=max_A, ncol=T_case) # Proportion susceptible
# exp_inc = matrix(NA,nrow=max_A, ncol=T_case) # Expected incidence
# exp_reported_cases = matrix(NA,nrow=AG, ncol=T_case) # Expected reported cases
# exp_inc_grouped = matrix(NA,nrow=AG, ncol=T_case) # Expected incidence in age groups
# 
# # build up the cumulative hazard for each individual based on their residence history
# # multiply by hr for the cluster and by hrs for covariates
# birthyear = data_prov$birthyear
# Y = data_prov$Y
# endres_history = as.matrix(serosurv_data[,paste0("endres_",2007:2023)])
# 
# beta_covs = 0
# 
# logit_lambda_end = -2
# logit_lambda_nonend = -2
# 
# 
# #Lambda by individual
# for (i in 1:data_prov$N) {
#   # Age 0 and 1 separately to account for possible maternal immunity
#   # Reduction factor
#   if (birthyear[i]>T_ss) { # for youngest children in serosurvey (without case data, we cannot estimate annual FOI in their first year of life)
#     lambda_i[i]=exp(log_cum_lambda)
#   } else if (birthyear[i]==T_ss) {
#     # for next youngest children in serosurvey (we can estimate FOI in their year of birth, but in no other years)
#     
#     # If resident of non-endemic country
#     if (endres_history[i,birthyear[i]]==0) {
#         lambda_i[i] = expit(logit_lambda_nonend)
#     }
# 
#     # If resident of endemic country
#     if (endres_history[i,birthyear[i]]==1) {
#       lambda_i[i] = expit(logit_lambda_end)
#     }
# 
#     # If resident of as
#     if (endres_history[i,birthyear[i]]==2) {
#       lambda_i[i] = lambda[birthyear[i]+max_A+T_case-T_ss]
#     }
# 
#     lambda_i[i] = lambda_i[i] + exp(log_cum_lambda) # add cumulative lambda from T_ss+1 to the year of the survey
# 
#   } else {
# 
#     # For all other children, we can estimate annual FOI from yob to T_ss
#     # Go through years 0 to T_ss
#     lambda_i[i] = 0
#     for (a in (birthyear[i]):T_ss) {
#       
#       # Note this will only add to individual FOI if endres_history[i,a] is 0, 1, or 2
#       # For years before child is born, endres_history is set to -1
# 
#       if (endres_history[i,a]==0) {
#         lambda_i[i] = lambda_i[i] + expit(logit_lambda_nonend)
#       }
# 
#       if (endres_history[i,a]==1) {
#         lambda_i[i] = lambda_i[i] + expit(logit_lambda_end)
#       }
# 
#       if (endres_history[i,a]==2) {
#         lambda_i[i] = lambda_i[i] + lambda[a+max_A+T_case-T_ss]
#       }
#     }
# 
#     if (T_ss<A_sp) {
#       lambda_i[i] = lambda_i[i] + exp(log_cum_lambda)
#     }
#   }
# 
#   # Expected probability of seropositivity given individual-level FOI (multiplied by HRs for covariates)
#   exp_sp[i] = 1-exp(-n_serotypes * lambda_i[i] * exp(sum(beta_covs*reg_vars[i,])))  
# }
# 
# # Now get hazard experienced by age cohorts for case data
# # First get population average lambda from serosurvey
# 
# # Sum of survey weights
# w = weights
# sumw=sum(w)
# 
# # Get lambda representative of AS population during serosurvey period by taking inverse-weighted sum
# # of FOIs from serosurvey participants
# for (t in 1:T_ss) {
# 
#   lambda_rep[t]=0
# 
#   for (i in 1:data_prov$N) {
#     # lambda_rep is representative FOI for AS based on serosurvey participants
#     lambda_rep[t] = lambda_rep[t] + lambda[t+max_A+T_case-T_ss] * exp(sum(beta_covs*reg_vars[i,])) * w[i] / sumw
#   }
# 
# }
# 
# # For the years not covered by the serosurvey, use lambda (from case data)
# for (t in 1:(max_A+T_case-T_ss)) {
#   lambda_temp[t] = lambda[t]
# }
# 
# # For the years covered by the serosurvey, use lambda_rep (from serosurvey participants)
# for (t in 1:T_ss) {
#   lambda_temp[t+max_A+T_case-T_ss] = lambda_rep[t]
# }
# 
# # Cumulative FOI by age in years (a) and year considered (t)
# cum_lambda[1, 1] = alpha * lambda_temp[max_A]
# for (a in 2:max_A) {
#   cum_lambda[a, 1] = sum(lambda_temp[(max_A-a+1):max_A])
# }
# for (t in 2:T_case) {
#   cum_lambda[1, t] = alpha * lambda_temp[max_A+t-1]#first age group born with no past exposure
#   for (a in 2:max_A) {
#     cum_lambda[a, t] = cum_lambda[a-1, t-1] + lambda_temp[max_A+t-1]
#   }
# }
# 
# # Expected case numbers (representing secondary cases)
# for (t in 1:T_case) {
#   for (a in 1:max_A) {
#     susc[a, t] = exp(-n_serotypes * cum_lambda[a, t]) 
#     mono[a, t] = n_serotypes * exp(-(n_serotypes-1) * cum_lambda[a, t]) * (1 - exp(-cum_lambda[a, t])) 
#     exp_inc[a, t] = (n_serotypes-1) * lambda_temp[max_A+t] * mono[a, t]
#   }
# }
# 
# # Expected reported case numbers by age groups in data (as expectation of negative binomial for the likelihood of observed case numbers)
# for (t in 1:T_case) {
#   for (a_gr in 1:AG) {
#     exp_inc_grouped[a_gr, t] = mean(exp_inc[lr_bound[a_gr]:ur_bound[a_gr], t])
#     exp_reported_cases[a_gr, t] = exp_inc_grouped[a_gr, t] * reporting_rate[t] * pop[a_gr, t]
#   }
# }

