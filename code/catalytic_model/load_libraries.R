require(dplyr)
require(rstan)
require(readxl)
require(tidyverse)
require(bayesplot)
require(binom)
require(data.table)
require(grid)
require(gridExtra)
require(ggpubr)
library(scales)  


# Function to estimate FOI based on mean age of infection
foi_init <- function(L, A, inc) {
      r0 <- 1 + L / A
      s <- 1 / r0
      inc / s
}

# Function to estimate incidence rate from case data
calc_incidence <- function(cases, pop, reporting_rate) {
      if (is.null(dim(pop))) {
            N <- sum(pop) # pop size
      } else {
            N <- sum(pop[, 1]) # pop size
      }
      total_cases <- sum(cases)
      nb_year <- ncol(cases) # years available of incidence data
      inc <- total_cases / nb_year / N / reporting_rate
      return(inc)
}


expit = function(x) {1/(1+exp(-x))}




run_FRAILTYregression_model_case_ss = function(compiled_model,
                                        serosurv_data,
                                        case_data,
                                        pop_data,
                                        lr_bound,
                                        ur_bound,
                                        svyweights,
                                        covariates,covariate_names,
                                        n_serotypes,runmod,suffix,
                                        firstyear_annualfoi,
                                        foi_sd = foi_sd,
                                        number_of_chains = 4,
                                        vary_reprateratio = 0,
                                        fixed_reprateratio = 1,
                                        combine_17_18=FALSE,
                                        lastcases=2023
                                        ) {
      
      # Define weights based on whether we are using survey weights or not
      serosurv_data$nonw = 1
      if (svyweights==1) {
            weights = 1/serosurv_data[,'prob']
            ## Normalize weights to sample size
            weights = weights / sum(weights) * nrow(serosurv_data)
      } else {
            weights = serosurv_data[,'nonw']
      }
      
      
      # Pick out covariates from serosurvey data based on names
      if (length(covariates)>0) {
            
            print(length(covariates))
            reg_vars=serosurv_data[,covariates,drop=FALSE]
      } else {
            reg_vars=array(0,dim=c(nrow(serosurv_data),0))
      }
      r_v = length(covariates)
      
      # Unique file suffix based on input parameters
      filesuffix = paste0('_fixeff_casess_covs_',paste0(covariates,collapse='.'),
                          '_ns_',n_serotypes,
                          suffix)
      
      print(filesuffix)
      
      
      
      #if 2017 and 2018 are combined into 1 year of transmission, change the yob
      earliest_yob <- 2007
      firstcases <-  2016
      if(combine_17_18==TRUE) {
            earliest_yob <- 2008
            firstcases <- 2017
            
      }
      
      # Input data for rstan
      data_prov <-
            list(
                  N = nrow(serosurv_data), # the number of individuals in the serosurvey
                  A_sp = 2023 - min(serosurv_data$yob) + 1, #number of age classes in the serosurvey
                  Y = serosurv_data$Seropositive, #serostatus of each individual
                  firstyear_exposure = serosurv_data$yob, # yob as calendar year
                  serosurv_year = 2023, # calendar year of the serosurvey
                  
                  w = unlist(weights),
                  endres_history = as.matrix(serosurv_data[,paste0("endres_",earliest_yob:2023)]),
                  r_v = r_v,
                  reg_vars = reg_vars,
                  
                  N_pos_control = N_pos_control,
                  control_tp = control_tp,
                  N_neg_control = N_neg_control,
                  control_fp = control_fp,
                  n_serotypes = n_serotypes,
                  
                  firstyear_annualfoi = firstyear_annualfoi,
                  
                  AG=AG,
                  max_A=max_A,
                  casedata_lastyear=lastcases,
                  casedata_firstyear=firstcases,
                  secondary_cases=case_data,
                  pop=pop_data,
                  lr_bound=lr_bound,
                  ur_bound=ur_bound,
                  foi_mean=foi_mean,
                  foi_sd_prior=foi_sd_prior,
                  foi_sd = foi_sd,
                  
                  reporting_mean=reporting_mean,
                  reporting_sd=reporting_sd,
                  
                  vary_reprateratio = vary_reprateratio,
                  fixed_reprateratio = fixed_reprateratio,
                  
                  old_serosurv_year = 2010,
                  old_serosurv_maxage = 40,
                  pos2010_18to25=179,
                  total2010_18to25 = 201,
                  pos2010_26to40 = 216,
                  total2010_26to40 = 217,
                  
                  clust=serosurv_data$cluster,
                  num_clust=length(unique(serosurv_data$cluster)),
                  alpha_priormean= 1.7,
                  alpha_priorvar= 1
                  
            )

      ### Without starting values
      if (runmod) {
            system.time(
                  as_denguefoi_fe <-
                        sampling(
                              compiled_model,
                              data = data_prov,
                              chains = number_of_chains,
                              iter = 2000,
                              warmup = 1000,
                              thin = 1,
                              control = list(adapt_delta = 0.80),
                              include = F,
                              pars = c('lambda_i','sumw','exp_sn','cum_lambda','exp_inc','exp_inc_grouped','exp_sp')
                        )
            )
            as_denguefoi_fe@stanmodel@dso <- new("cxxdso")
            saveRDS(as_denguefoi_fe, file = paste0('as_denguefoi',filesuffix,'.rds'))
            
      } else {
            as_denguefoi_fe = readRDS(paste0('as_denguefoi',filesuffix,'.rds'))
      }

      return(as_denguefoi_fe)
      
}




### graph 1: seroprevalence fit graph

seroprev_fitgraph <- function(serosurv_data, as_denguefoi_fe){
      
      
      exp_obs_sp = rstan::extract(as_denguefoi_fe,pars="exp_obs_sp", inc_warmup=FALSE, permute=T)$exp_obs_sp
      exp_obs_sp_mean <- apply(exp_obs_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
      
      exp_obs_sp_byage = sapply(1:2000,function(i) {
            serosurv_data$exp_obs_sptemp=exp_obs_sp[i,]
            exp_obs_sp_byage = serosurv_data %>% group_by(yob) %>% summarise(esp=mean(exp_obs_sptemp))
            return(exp_obs_sp_byage$esp)
      })
      exp_obs_sp_byage_summ = apply(exp_obs_sp_byage,1,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      lambda = rstan::extract(as_denguefoi_fe,pars="lambda_temp",inc_warmup=F,permute=T)$lambda_temp
      lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
      colnames(lambda)=1967 + (0:(ncol(lambda)-1))
      
      pre_lambda = rstan::extract(as_denguefoi_fe,pars="pre_lambda",inc_warmup=F,permute=F)
      pre_lambda = apply(pre_lambda,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      post_lambda_sp = rstan::extract(as_denguefoi_fe,pars="post_lambda",inc_warmup=F,permute=F)
      post_lambda_sp = apply(post_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
      end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
      nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      bb_data = serosurv_data %>% group_by(yob) %>% summarise(n=n(),x=sum(Seropositive))
      bin_cis = as.data.frame(t(sapply(1:nrow(bb_data),
                                       function(x) {
                                             d=binom.confint(bb_data$x[x],bb_data$n[x])
                                             return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
                                       })))
      colnames(bin_cis)=c("SP","lwr","upr")
      bin_cis$yob = bb_data$yob
      
      d_forplot = data.frame('Age'=rep(bin_cis$yob,2),
                             'Model'=rep(c('Observed','Model estimated'),each=nrow(bin_cis)),
                             'SP'=c(bin_cis[,1],
                                    exp_obs_sp_byage_summ[2,]),
                             'lwr'=c(bin_cis[,2],exp_obs_sp_byage_summ[1,]),
                             'upr'=c(bin_cis[,3],exp_obs_sp_byage_summ[3,]))
      
      theme_set(theme_classic())
      psero=ggplot() +
            geom_errorbar(data = d_forplot %>% filter(Model=="Observed"),
                          aes(x=2023-Age+1,ymin=lwr,ymax=upr),col='black',width=0) +
            geom_line(data = d_forplot,aes(x=2023-Age+1,y=SP,col=factor(Model,levels=c('Observed','Model estimated'))))+
            geom_ribbon(data = d_forplot,aes(x=2023-Age+1,ymin=lwr,ymax=upr,fill=factor(Model,levels=c('Observed','Model estimated'))),alpha=0.3)+
            geom_point(data = d_forplot, aes(x=2023-Age+1,y=SP,shape=factor(Model,levels=c('Observed','Model estimated'))))+
            scale_color_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
            scale_fill_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
            scale_shape_manual(name=NULL,values=c(16,NA),labels=c('Observed','Model estimated'))+
            ylab('Seroprevalence')+ylim(c(0,1))+
            scale_x_continuous(name='Birth year',breaks=9:17,labels=2015:2007)+
            theme(axis.text=element_text(size=6),
                  axis.title = element_text(size=6),
                  legend.title = element_text(size=6),
                  legend.text = element_text(size=6),
                  legend.position = c(0.1,0.89),
                  strip.text = element_text(size=6))
      
      return(psero)
      
}


### graph 2: year fit graph

year_fitgraph <- function(as_denguefoi_fe,yr1=2016,lastyr=2023){
      

      reprate = rstan::extract(as_denguefoi_fe,'reporting_rate',inc_warmup=F,permute=T)$reporting_rate
      reprate = quantile(reprate,c(0.025,0.5,0.975))
      
      susc = rstan::extract(as_denguefoi_fe,'susc',inc_warmup=F,permute=T)$susc
      susc = apply(susc,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
      
      mono = rstan::extract(as_denguefoi_fe,'mono',inc_warmup=F,permute=T)$mono
      mono = apply(mono,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
      
      repcases = rstan::extract(as_denguefoi_fe,'exp_reported_cases')$exp_reported_cases
      repcases = apply(repcases,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
      repcases_forplot = as.data.frame(apply(repcases,1,cbind))
      colnames(repcases_forplot)=c("LCI","Median","UCI")
      repcases_forplot$AgeGroup=agegroups
      repcases_forplot$Year=rep(yr1:lastyr,each=10)
      
      repcasestotal_byyear = apply(repcases,c(1,3),sum)
      repcasestotal_byyear = as.data.frame(t(apply(repcasestotal_byyear,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
      colnames(repcasestotal_byyear)=c("LCI","Median","UCI")
      repcasestotal_byyear$Year=yr1:lastyr
      repcasestotal_byyear$Model='Estimated'
      
      obscasestotal = data.frame('Median'=colSums(case_data),
                                 'LCI'=NA,
                                 'UCI'=NA,
                                 'Year'=yr1:lastyr,
                                 'Model'="Observed")
      repcasestotal_byyear = rbind(repcasestotal_byyear,obscasestotal)
      
      
      p1=ggplot() + geom_col(data=repcasestotal_byyear %>% filter(Model=="Observed"),aes(x=Year,y=Median),col='grey',alpha=.3)+
            geom_point(data=repcasestotal_byyear %>% filter(Model=="Estimated"),aes(x = Year, y = Median),col = "steelblue") +
            geom_errorbar(data=repcasestotal_byyear %>% filter(Model=="Estimated"),
                          aes(x=Year,ymin=LCI,ymax=UCI),col = "steelblue")+
            scale_x_continuous(name="Year",breaks=seq(2010,lastyr,by=2))+ylab("Reported cases")+
            theme(axis.text=element_text(size=6),
                  axis.title=element_text(size=6))
      
      
      return(p1)
      
}



### graph 3: age group fit graph

agegroup_fitgraph <- function(as_denguefoi_fe = as_denguefoi_fe){
      
      repcases = rstan::extract(as_denguefoi_fe,'exp_reported_cases')$exp_reported_cases
      repcases = apply(repcases,c(2,3),function(x) quantile(x,c(0.025,0.5,0.975)))
      
      repcasestotal_byage = apply(repcases,c(1,2),sum)
      repcasestotal_byage = as.data.frame(t(apply(repcasestotal_byage,2,function(x) quantile(x,c(0.025,0.5,0.975)))))
      colnames(repcasestotal_byage)=c("LCI","Median","UCI")
      repcasestotal_byage$AgeGroup=agegroups
      repcasestotal_byage$Model='Estimated'
      
      obscasestotal = data.frame('Median'=rowSums(case_data),
                                 'LCI'=NA,
                                 'UCI'=NA,
                                 'AgeGroup'=agegroups,
                                 'Model'="Observed")
      repcasestotal_byage = rbind(repcasestotal_byage,obscasestotal)
      
      repcasestotal_byage$AgeGroup = factor(repcasestotal_byage$AgeGroup,
                                            levels=agegroups)
      
      agegroups_plot = c("1-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49")
      p2=ggplot() + 
            geom_col(data=repcasestotal_byage %>% 
                                   filter(Model=="Observed"),
                             aes(x=AgeGroup,y=Median),col='grey',alpha=.3)+
            geom_point(data=repcasestotal_byage %>% filter(Model=="Estimated"),aes(x = AgeGroup, y = Median),col = "steelblue") +
            geom_errorbar(data=repcasestotal_byage %>% filter(Model=="Estimated"),
                          aes(x=AgeGroup,ymin=LCI,ymax=UCI),col = "steelblue")+
            scale_x_discrete(name="Age Group",labels=agegroups_plot,guide = guide_axis(n.dodge = 2))+
            ylab("Reported cases")+
            theme(axis.text=element_text(size=6),
                  axis.title=element_text(size=6))
      
      p2
      
      
}


make_foi_table <- function(as_denguefoi_fe,lastyear_annualfoi,firstyear_annualfoi,serosurv_year,firstyear_toreport){
      
      lambda = rstan::extract(as_denguefoi_fe,pars="lambda_temp",inc_warmup=F,permute=T)$lambda_temp
      lambdapos = apply(lambda,2,function(x) mean(x>1e-3))
      names(lambdapos) = 1967:serosurv_year
      lambdapos = lambdapos[names(lambdapos) %in% firstyear_annualfoi:serosurv_year]
      #lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
      lambda = apply(lambda,2,function(x) c(quantile(x,c(0.025)),mean(x),quantile(x,0.975)))
      
      if (lastyear_annualfoi<serosurv_year) {
        post_lambda = rstan::extract(as_denguefoi_fe,pars="post_lambda",inc_warmup=F,permute=F)
        post_lambda = apply(post_lambda,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      }
      
      pre_lambda = rstan::extract(as_denguefoi_fe,pars="pre_lambda",inc_warmup=F,permute=F)
      pre_lambda = apply(pre_lambda,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
      end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
      nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
      
      if (lastyear_annualfoi<serosurv_year) {
        fois <- c(     ###American samoa
              paste0(round(n_serotypes*lambda[2,],2),' (',
                     round(n_serotypes*lambda[1,],2),',',
                     round(n_serotypes*lambda[3,],2),')'),
              #### from the last year we can estiamte annual FOI to the year of the serosurvey
              rep(paste0(round(n_serotypes*post_lambda[2],2),' (',
                         round(n_serotypes*post_lambda[1],2),',',
                         round(n_serotypes*post_lambda[3],2),')'),1),
              #### from before the first year we can estimate annual FOI 
              rep(paste0(round(n_serotypes*pre_lambda[2],2),' (',
                         round(n_serotypes*pre_lambda[1],2),',',
                         round(n_serotypes*pre_lambda[3],2),')'),1),
              ### endemic
              paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
                     round(n_serotypes*end_lambda_sp[1],2),',',
                     round(n_serotypes*end_lambda_sp[3],2),')'),
              
              ### non endemic
              paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
                     round(n_serotypes*nonend_lambda_sp[1],2),',',
                     round(n_serotypes*nonend_lambda_sp[3],2),')')
              
        )
        name=c((lastyear_annualfoi-length(fois)+3):lastyear_annualfoi,
               paste0("Annual FOI ",lastyear_annualfoi+1,"-",serosurv_year),
               paste0("Annual FOI 1967-",firstyear_annualfoi),
               "Endemic","Non-Endemic")
      } else {
        fois <- c(     ###American samoa
          paste0(round(n_serotypes*lambda[2,],2),' (',
                 round(n_serotypes*lambda[1,],2),',',
                 round(n_serotypes*lambda[3,],2),')'),
          #### from before the first year we can estimate annual FOI 
          rep(paste0(round(n_serotypes*pre_lambda[2],2),' (',
                     round(n_serotypes*pre_lambda[1],2),',',
                     round(n_serotypes*pre_lambda[3],2),')'),1),
          ### endemic
          paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
                 round(n_serotypes*end_lambda_sp[1],2),',',
                 round(n_serotypes*end_lambda_sp[3],2),')'),
          
          ### non endemic
          paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
                 round(n_serotypes*nonend_lambda_sp[1],2),',',
                 round(n_serotypes*nonend_lambda_sp[3],2),')')
        )
        name=c((lastyear_annualfoi-length(fois)+4):lastyear_annualfoi,
               paste0("Annual FOI 1967-",firstyear_annualfoi),
               "Endemic","Non-Endemic")
      }
      
      length(fois)
      
      
      
      out_table <- data.frame(name,fois) %>%
            filter(!name %in%(lastyear_annualfoi-length(fois)+3):(firstyear_toreport-1))
      return(out_table)
      
      
}

seroprev_fitgraph_byschool <- function(serosurv_data, as_denguefoi_fe){
  
  
  exp_obs_sp = rstan::extract(as_denguefoi_fe,pars="exp_obs_sp", inc_warmup=FALSE, permute=T)$exp_obs_sp
  exp_obs_sp_mean <- apply(exp_obs_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
  
  exp_obs_sp_byageschool = sapply(1:2000,function(i) {
    serosurv_data$exp_obs_sptemp=exp_obs_sp[i,]
    exp_obs_sp_byage = serosurv_data %>% group_by(yob,School) %>% summarise(esp=mean(exp_obs_sptemp),.groups="drop_last")
    return(exp_obs_sp_byage$esp)
  })
  exp_obs_sp_byageschool_summ = apply(exp_obs_sp_byageschool,1,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  bb_data = serosurv_data %>% group_by(yob,School) %>% summarise(n=n(),x=sum(Seropositive),.groups="drop_last")
  bin_cis = as.data.frame(t(sapply(1:nrow(bb_data),
                                   function(x) {
                                     d=binom.confint(bb_data$x[x],bb_data$n[x])
                                     return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
                                   })))
  colnames(bin_cis)=c("SP","lwr","upr")
  bin_cis$yob = bb_data$yob
  
  d_forplot = data.frame('Age'=rep(bin_cis$yob,2),
                         'School'=rep(bb_data$School,2),
                         'Model'=rep(c('Observed','Model estimated'),each=nrow(bin_cis)),
                         'SP'=c(bin_cis[,1],
                                exp_obs_sp_byageschool_summ[2,]),
                         'lwr'=c(bin_cis[,2],exp_obs_sp_byageschool_summ[1,]),
                         'upr'=c(bin_cis[,3],exp_obs_sp_byageschool_summ[3,]))
  
  theme_set(theme_classic())
  psero=ggplot() +
    geom_errorbar(data = d_forplot %>% filter(Model=="Observed"),
                  aes(x=2023-Age+1,ymin=lwr,ymax=upr),col='black',width=0) +
    geom_line(data = d_forplot,aes(x=2023-Age+1,y=SP,col=factor(Model,levels=c('Observed','Model estimated'))))+
    geom_ribbon(data = d_forplot,aes(x=2023-Age+1,ymin=lwr,ymax=upr,fill=factor(Model,levels=c('Observed','Model estimated'))),alpha=0.3)+
    geom_point(data = d_forplot, aes(x=2023-Age+1,y=SP,shape=factor(Model,levels=c('Observed','Model estimated'))))+
    scale_color_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_fill_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_shape_manual(name=NULL,values=c(16,NA),labels=c('Observed','Model estimated'))+
    ylab('Seroprevalence')+ylim(c(0,1))+
    scale_x_continuous(name='Birth year',breaks=9:17,labels=2015:2007)+
    facet_wrap(vars(School)) +
    theme(axis.text=element_text(size=6),
          axis.title = element_text(size=6),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position = c(0.1,0.89),
          strip.text = element_text(size=6))
  
  return(psero)
  
}



### Serosurvey only
run_FRAILTYregression_model_ss = function(compiled_model,
                                          serosurv_data,
                                          svyweights,
                                          covariates,covariate_names,
                                          n_serotypes,runmod,suffix,
                                          foi_sd = foi_sd,
                                          number_of_chains = 4,
                                          combine_17_18=FALSE
) {
  
  # Define weights based on whether we are using survey weights or not
  serosurv_data$nonw = 1
  if (svyweights==1) {
    weights = 1/serosurv_data[,'prob']
    ## Normalize weights to sample size
    weights = weights / sum(weights) * nrow(serosurv_data)
  } else {
    weights = serosurv_data[,'nonw']
  }
  
  
  # Pick out covariates from serosurvey data based on names
  if (length(covariates)>0) {
    
    print(length(covariates))
    reg_vars=serosurv_data[,covariates,drop=FALSE]
  } else {
    reg_vars=array(0,dim=c(nrow(serosurv_data),0))
  }
  r_v = length(covariates)
  
  # Unique file suffix based on input parameters
  filesuffix = paste0('_fixeff_ssonly_covs_',paste0(covariates,collapse='.'),
                      '_ns_',n_serotypes,
                      suffix)
  
  print(filesuffix)
  
  
  
  #if 2017 and 2018 are combined into 1 year of transmission, change the yob
  earliest_yob <- 2007
  if(combine_17_18==TRUE) {
    earliest_yob <- 2008
  }
  
  # Input data for rstan
  data_prov <-
    list(
      N = nrow(serosurv_data), # the number of individuals in the serosurvey
      A_sp = 2023 - min(serosurv_data$yob) + 1, #number of age classes in the serosurvey
      Y = serosurv_data$Seropositive, #serostatus of each individual
      firstyear_exposure = serosurv_data$yob, # yob as calendar year
      serosurv_year = 2023, # calendar year of the serosurvey
      
      w = unlist(weights),
      endres_history = as.matrix(serosurv_data[,paste0("endres_",earliest_yob:2023)]),
      r_v = r_v,
      reg_vars = reg_vars,
      
      N_pos_control = N_pos_control,
      control_tp = control_tp,
      N_neg_control = N_neg_control,
      control_fp = control_fp,
      n_serotypes = n_serotypes,
      
      foi_mean=foi_mean,
      foi_sd_prior=foi_sd_prior,
      foi_sd = foi_sd,
      
      clust=serosurv_data$cluster,
      num_clust=length(unique(serosurv_data$cluster)),
      alpha_priormean= 1.7,
      alpha_priorvar= 1
      
    )
  
  ### Without starting values
  if (runmod) {
    system.time(
      as_denguefoi_fe <-
        sampling(
          compiled_model,
          data = data_prov,
          chains = number_of_chains,
          iter = 2000,
          warmup = 1000,
          thin = 1,
          control = list(adapt_delta = 0.80),
          include = F,
          pars = c('lambda_i','sumw','exp_sn','cum_lambda','exp_inc','exp_inc_grouped','exp_sp')
        )
    )
    as_denguefoi_fe@stanmodel@dso <- new("cxxdso")
    saveRDS(as_denguefoi_fe, file = paste0('as_denguefoi',filesuffix,'.rds'))
    
  } else {
    as_denguefoi_fe = readRDS(paste0('as_denguefoi',filesuffix,'.rds'))
  }
  
  return(as_denguefoi_fe)
  
}

### graph 1: seroprevalence fit graph

seroprev_fitgraph_ssonly <- function(serosurv_data, as_denguefoi_fe){
  
  
  exp_obs_sp = rstan::extract(as_denguefoi_fe,pars="exp_obs_sp", inc_warmup=FALSE, permute=T)$exp_obs_sp
  exp_obs_sp_mean <- apply(exp_obs_sp, 2, function(x) quantile(x,c(0.025,0.5,0.975))) # mean across chains and run
  
  exp_obs_sp_byage = sapply(1:2000,function(i) {
    serosurv_data$exp_obs_sptemp=exp_obs_sp[i,]
    exp_obs_sp_byage = serosurv_data %>% group_by(yob) %>% summarise(esp=mean(exp_obs_sptemp))
    return(exp_obs_sp_byage$esp)
  })
  exp_obs_sp_byage_summ = apply(exp_obs_sp_byage,1,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  lambda = rstan::extract(as_denguefoi_fe,pars="lambda",inc_warmup=F,permute=T)$lambda
  lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
  colnames(lambda)=2007 + (0:(ncol(lambda)-1))
  
  post_lambda_sp = rstan::extract(as_denguefoi_fe,pars="post_lambda",inc_warmup=F,permute=F)
  post_lambda_sp = apply(post_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
  end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
  nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  bb_data = serosurv_data %>% group_by(yob) %>% summarise(n=n(),x=sum(Seropositive))
  bin_cis = as.data.frame(t(sapply(1:nrow(bb_data),
                                   function(x) {
                                     d=binom.confint(bb_data$x[x],bb_data$n[x])
                                     return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
                                   })))
  colnames(bin_cis)=c("SP","lwr","upr")
  bin_cis$yob = bb_data$yob
  
  d_forplot = data.frame('Age'=rep(bin_cis$yob,2),
                         'Model'=rep(c('Observed','Model estimated'),each=nrow(bin_cis)),
                         'SP'=c(bin_cis[,1],
                                exp_obs_sp_byage_summ[2,]),
                         'lwr'=c(bin_cis[,2],exp_obs_sp_byage_summ[1,]),
                         'upr'=c(bin_cis[,3],exp_obs_sp_byage_summ[3,]))
  
  theme_set(theme_classic())
  psero=ggplot() +
    geom_errorbar(data = d_forplot %>% filter(Model=="Observed"),
                  aes(x=2023-Age+1,ymin=lwr,ymax=upr),col='black',width=0) +
    geom_line(data = d_forplot,aes(x=2023-Age+1,y=SP,col=factor(Model,levels=c('Observed','Model estimated'))))+
    geom_ribbon(data = d_forplot,aes(x=2023-Age+1,ymin=lwr,ymax=upr,fill=factor(Model,levels=c('Observed','Model estimated'))),alpha=0.3)+
    geom_point(data = d_forplot, aes(x=2023-Age+1,y=SP,shape=factor(Model,levels=c('Observed','Model estimated'))))+
    scale_color_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_fill_manual(name=NULL,values=c(NA,'blue'),na.value=NA,labels=c('Observed','Model estimated'))+
    scale_shape_manual(name=NULL,values=c(16,NA),labels=c('Observed','Model estimated'))+
    ylab('Seroprevalence')+ylim(c(0,1))+
    scale_x_continuous(name='Birth year',breaks=9:17,labels=2015:2007)+
    theme(axis.text=element_text(size=6),
          axis.title = element_text(size=6),
          legend.title = element_text(size=6),
          legend.text = element_text(size=6),
          legend.position = c(0.1,0.89),
          strip.text = element_text(size=6))
  
  return(psero)
  
}

make_foi_table_ssonly <- function(as_denguefoi_fe,serosurv_year){
  
  lambda = rstan::extract(as_denguefoi_fe,pars="lambda",inc_warmup=F,permute=T)$lambda
  lambda = apply(lambda,2,function(x) c(quantile(x,c(0.025)),mean(x),quantile(x,0.975)))
  
  if (lastyear_annualfoi<serosurv_year) {
    post_lambda = rstan::extract(as_denguefoi_fe,pars="post_lambda",inc_warmup=F,permute=F)
    post_lambda = apply(post_lambda,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  }
  
  end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
  end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
  nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  fois <- c(     ###American samoa
    paste0(round(n_serotypes*lambda[2,],2),' (',
           round(n_serotypes*lambda[1,],2),',',
           round(n_serotypes*lambda[3,],2),')'),
    #### from the last year we can estiamte annual FOI to the year of the serosurvey
    rep(paste0(round(n_serotypes*post_lambda[2],2),' (',
               round(n_serotypes*post_lambda[1],2),',',
               round(n_serotypes*post_lambda[3],2),')'),1),
    ### endemic
    paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
           round(n_serotypes*end_lambda_sp[1],2),',',
           round(n_serotypes*end_lambda_sp[3],2),')'),
    
    ### non endemic
    paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
           round(n_serotypes*nonend_lambda_sp[1],2),',',
           round(n_serotypes*nonend_lambda_sp[3],2),')')
    
  )
  name=c(serosurv_firstyear:(serosurv_lastyear-1),
         paste0("Annual FOI ",serosurv_lastyear,"-",serosurv_year),
         "Endemic","Non-Endemic")
  
  out_table <- data.frame(name,fois)
  
  return(out_table)
  
  
}

run_FRAILTYregression_model_case_ss_old = function(compiled_model,
                                               serosurv_data,
                                               case_data,
                                               pop_data,
                                               lr_bound,
                                               ur_bound,
                                               svyweights,
                                               covariates,covariate_names,
                                               n_serotypes,runmod,suffix,
                                               number_of_chains = 4,
                                               combine_17_18=FALSE,
                                               lastcases=2023
) {
  
  # Define weights based on whether we are using survey weights or not
  serosurv_data$nonw = 1
  if (svyweights==1) {
    weights = 1/serosurv_data[,'prob']
    ## Normalize weights to sample size
    weights = weights / sum(weights) * nrow(serosurv_data)
  } else {
    weights = serosurv_data[,'nonw']
  }
  
  
  # Pick out covariates from serosurvey data based on names
  if (length(covariates)>0) {
    
    print(length(covariates))
    reg_vars=serosurv_data[,covariates,drop=FALSE]
  } else {
    reg_vars=array(0,dim=c(nrow(serosurv_data),0))
  }
  r_v = length(covariates)
  
  # Unique file suffix based on input parameters
  filesuffix = paste0('_fixeff_casess_covs_',paste0(covariates,collapse='.'),
                      '_ns_',n_serotypes,
                      suffix)
  
  print(filesuffix)
  
  
  
  #if 2017 and 2018 are combined into 1 year of transmission, change the yob
  earliest_yob <- 2007
  firstcases <-  2016
  if(combine_17_18==TRUE) {
    earliest_yob <- 2008
    firstcases <- 2017
    
  }
  
  # Input data for rstan
  data_prov <-
    list(
      N = nrow(serosurv_data), # the number of individuals in the serosurvey
      A_sp = 2023 - min(serosurv_data$yob) + 1, #number of age classes in the serosurvey
      Y = serosurv_data$Seropositive, #serostatus of each individual
      w = unlist(weights),
      
      endres_history = as.matrix(serosurv_data[,paste0("endres_",earliest_yob:2023)]),
      r_v = r_v,
      reg_vars = reg_vars,
      
      birthyear = serosurv_data$yob, # yob as calendar year
      birthmonth = rep(1,nrow(serosurv_data)),
      serosurv_year = 2023, # calendar year of the serosurvey
      
      
      N_pos_control = N_pos_control,
      control_tp = control_tp,
      N_neg_control = N_neg_control,
      control_fp = control_fp,
      n_serotypes = n_serotypes,
      
      AG=AG,
      max_A=max_A,
      casedata_lastyear=lastcases,
      casedata_firstyear=firstcases,
      secondary_cases=case_data,
      pop=pop_data,
      lr_bound=lr_bound,
      ur_bound=ur_bound,
      foi_mean=foi_mean,
      foi_sd = foi_sd,
      
      reporting_mean=reporting_mean,
      reporting_sd=reporting_sd,
      old_serosurv_year = 2010,
      old_serosurv_maxage = 40,
      pos2010_18to25=179,
      total2010_18to25 = 201,
      pos2010_26to40 = 216,
      total2010_26to40 = 217,
      
      clust=serosurv_data$cluster,
      num_clust=length(unique(serosurv_data$cluster)),
      alpha_priormean= 1.7,
      alpha_priorvar= 1
      
    )
  
  ### Without starting values
  if (runmod) {
    system.time(
      as_denguefoi_fe <-
        sampling(
          compiled_model,
          data = data_prov,
          chains = number_of_chains,
          iter = 1000,
          warmup = 500,
          thin = 1,
          control = list(adapt_delta = 0.80),
          include = F,
          pars = c('lambda_i','sumw','exp_sn','cum_lambda','exp_inc','exp_inc_grouped','exp_sp')
        )
    )
    as_denguefoi_fe@stanmodel@dso <- new("cxxdso")
    saveRDS(as_denguefoi_fe, file = paste0('as_denguefoi',filesuffix,'.rds'))
    
  } else {
    as_denguefoi_fe = readRDS(paste0('as_denguefoi',filesuffix,'.rds'))
  }
  
  return(as_denguefoi_fe)
  
}


make_foi_table_old <- function(as_denguefoi_fe,lastyear_annualfoi,serosurv_year,firstyear_toreport){
  
  lambda = rstan::extract(as_denguefoi_fe,pars="lambda_temp",inc_warmup=F,permute=T)$lambda_temp
  lambdapos = apply(lambda,2,function(x) mean(x>1e-3))
  names(lambdapos) = 1967:serosurv_year
  lambdapos = lambdapos[names(lambdapos) %in% firstyear_toreport:serosurv_year]
  #lambda = apply(lambda,2,function(x) quantile(x,c(0.025,0.5,0.975)))
  lambda = apply(lambda,2,function(x) c(quantile(x,c(0.025)),mean(x),quantile(x,0.975)))
  
  if (lastyear_annualfoi<serosurv_year) {
    post_lambda = rstan::extract(as_denguefoi_fe,pars="post_lambda",inc_warmup=F,permute=F)
    post_lambda = apply(post_lambda,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  }
  
  end_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_end",inc_warmup=F,permute=F))
  end_lambda_sp = apply(end_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  nonend_lambda_sp = expit(rstan::extract(as_denguefoi_fe,pars="logit_lambda_nonend",inc_warmup=F,permute=F))
  nonend_lambda_sp = apply(nonend_lambda_sp,3,function(x) quantile(x,c(0.025,0.5,0.975)))
  
  if (lastyear_annualfoi<serosurv_year) {
    fois <- c(     ###American samoa
      paste0(round(n_serotypes*lambda[2,],2),' (',
             round(n_serotypes*lambda[1,],2),',',
             round(n_serotypes*lambda[3,],2),')'),
      #### from the last year we can estiamte annual FOI to the year of the serosurvey
      rep(paste0(round(n_serotypes*post_lambda[2],2),' (',
                 round(n_serotypes*post_lambda[1],2),',',
                 round(n_serotypes*post_lambda[3],2),')'),1),
      ### endemic
      paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
             round(n_serotypes*end_lambda_sp[1],2),',',
             round(n_serotypes*end_lambda_sp[3],2),')'),
      
      ### non endemic
      paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
             round(n_serotypes*nonend_lambda_sp[1],2),',',
             round(n_serotypes*nonend_lambda_sp[3],2),')')
      
    )
    name=c((lastyear_annualfoi-length(fois)+3):lastyear_annualfoi,
           paste0("Annual FOI ",lastyear_annualfoi+1,"-",serosurv_year),
           "Endemic","Non-Endemic")
  } else {
    fois <- c(     ###American samoa
      paste0(round(n_serotypes*lambda[2,],2),' (',
             round(n_serotypes*lambda[1,],2),',',
             round(n_serotypes*lambda[3,],2),')'),
      ### endemic
      paste0(round(n_serotypes*end_lambda_sp[2],2),' (',
             round(n_serotypes*end_lambda_sp[1],2),',',
             round(n_serotypes*end_lambda_sp[3],2),')'),
      
      ### non endemic
      paste0(round(n_serotypes*nonend_lambda_sp[2],2),' (',
             round(n_serotypes*nonend_lambda_sp[1],2),',',
             round(n_serotypes*nonend_lambda_sp[3],2),')')
    )
    name=c((lastyear_annualfoi-length(fois)+3):lastyear_annualfoi,
           "Endemic","Non-Endemic")
  }
  
  length(fois)
  
  
  
  out_table <- data.frame(name,fois) %>%
    filter(!name %in%(lastyear_annualfoi-length(fois)+3):(firstyear_toreport-1))
  return(out_table)
  
  
}

