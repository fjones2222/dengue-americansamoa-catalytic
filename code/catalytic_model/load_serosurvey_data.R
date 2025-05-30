

survweights <- read_csv('Data/AS Data_Probability.csv')


wide_dataset <- read_csv("Data/wide_dataset.csv")%>% 
      rename(RecordID=`Record ID`,
               Age=`Calculated age based on DOB.`,
               Village=`Village of Residence`,
             Date.of.birth=`Date of birth`
             ) %>%
      mutate(Seropositive = as.numeric(`Final Test Result`=="P"),
             sexf = as.numeric(Sex=="Female"),
             cluster=as.numeric(as.factor(School))
             ) %>%
      left_join(survweights,by='RecordID')
     

long_dataset <- read_csv("Data/long_dataset.csv")


analysis_set <- list()

### dataset 1: everyone lived in AS all the time

analysis_set[["everyone"]] <- long_dataset %>%
      select(RecordID=`Record ID`, Latter_location_group,year) %>%
      mutate(Latter_location_group="American Samoa") %>%
      arrange(RecordID,year)%>%
      mutate(year=paste0("res_",year))%>%
      spread(year,Latter_location_group) %>%
      right_join(wide_dataset)


### dataset 2: only use data from people who lived in AS all the time

analysis_set[["lifelong"]]<- long_dataset %>%
      select(RecordID=`Record ID`, Latter_location_group,year) %>%
      arrange(RecordID,year)%>%
      mutate(year=paste0("res_",year))%>%
      spread(year,Latter_location_group) %>%
      right_join(wide_dataset) #%>%
      #filter(`5. My child has been a resident of American Samoa`=="All his/her life")

analysis_set[["lifelong"]]$lifelong = 
  as.numeric(apply(analysis_set[["lifelong"]],1,function(x) sum(!is.na(x[paste0('res_',2007:2023)]) & x[paste0('res_',2007:2023)] != "American Samoa")==0))

### dataset 3: use the former location

analysis_set[["former"]] <-long_dataset %>% select(RecordID=`Record ID`, Former_location_group,year) %>%
      arrange(RecordID,year)%>%
      mutate(year=paste0("res_",year))%>%
      spread(year,Former_location_group) %>%
      right_join(wide_dataset)

### dataset 4: use the latter location

analysis_set[["latter"]]<-long_dataset %>% select(RecordID=`Record ID`, Latter_location_group,year) %>%
      arrange(RecordID,year)%>%
      mutate(year=paste0("res_",year))%>%
      spread(year,Latter_location_group) %>%
      right_join(wide_dataset)



#make the categories proper for each of the correct numbers
for(name in names(analysis_set)){
      

      analysis_set[[name]] <- analysis_set[[name]] %>% 
            mutate(across(res_2007:res_2023,.names="end{.col}",
                          .fns=~case_when(
                                is.na(.x) ~ -1,
                                .x=="Non-endemic" ~ 0,
                                .x=="Pacific Islands" ~ 1,
                                .x=="American Samoa" ~ 2
                                
                          ))) 
                         
      
      
}


# analysis_set$everyone$endres_2007 %>% table()




# 
# ## Serosurvey data
# serosurv_data = read.csv('Data/merged_latter1.csv') %>% dplyr::select(-X) %>% rename(RecordID=Record.ID,
#                                                                                      Age=Calculated.age.based.on.DOB.,
#                                                                                      Village=Village.of.Residence) %>%
#       filter(Final.Test.Result!="I")#%>%
#       # filter(!RecordID %in% c(799, 197, 174,850, 792, 415, 251, 187, 825, 778, 296 ))
#       
# colnames(serosurv_data) = gsub('X','res_',colnames(serosurv_data))
# 
# 
# # serosurv_data = 
# # serosurv_data = serosurv_data 
# ## Define endemicity of residence by year
# allcountries = unique(unlist(c(serosurv_data[,paste0('res_',2007:2022)])))
# endemic_list = c("Fiji","Philippines","Western Samoa","Tuvalu")
# 
# for (y in 2007:2023) {
#       
#       serosurv_data[,paste0("endres_",y)] = -1*as.numeric(is.na(serosurv_data[,paste0("res_",y)])) +
#             as.numeric(serosurv_data[,paste0("res_",y)] %in% endemic_list) +
#             2*as.numeric(serosurv_data[,paste0("res_",y)] %in% c("American Samoa"))
#       
#       
#       # serosurv_data[,paste0("endres_",y)] <- -1*as.numeric(is.na(serosurv_data[,paste0("res_",y)]))
#       
#       # #everyone lives in AS
#       # serosurv_data[,paste0("endres_",y)] = -1*as.numeric(is.na(serosurv_data[,paste0("res_",y)])) + 
#       #       # as.numeric(serosurv_data[,paste0("res_",y)] %in% endemic_list) + 
#       #       2*as.numeric(!is.na(serosurv_data[,paste0("res_",y)]))
#       
# }





# # Observed seroprevalence and binomial CIs by birth year and island
# bb_data = serosurv_data %>% group_by(yob) %>% summarise(n=n(),x=sum(Seropositive))
# bin_cis = t(sapply(1:nrow(bb_data),
#                    function(x) {
#                          d=binom.confint(bb_data$x[x],bb_data$n[x])
#                          return(rbind(d$mean[d$method=="wilson"],d$lower[d$method=="wilson"],d$upper[d$method=="wilson"]))
#                    }))
# colnames(bin_cis)=c("SP","lwr","upr")




