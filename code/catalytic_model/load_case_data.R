# Case data
case_data_indiv = xl.read.file("data/raw_data/ArboNET Case Data.xlsx")
case_data_indiv$age = as.numeric(case_data_indiv$age)

cases_by_year = case_data_indiv %>% group_by(year) %>% summarise(ncase = n())
cases_by_year = rbind(cases_by_year,data.frame('year'=c(2019:2021,2023),'ncase'=0)) %>% arrange(year)
# old cases
# cases_by_year = rbind(cases_by_year,data.frame('year'=c(2010:2015,2019:2021,2023),'ncase'=0)) %>% arrange(year)



allagegroups = paste0(seq(0,80,by=5)," to ",seq(4,84,by=5))
agegroups = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")

caseage_by_year = case_data_indiv %>% filter(age<50 & age>0) %>%
      mutate(AgeGroup = cut(age,breaks=seq(-0.01,50,by=5),labels=agegroups)) %>% 
      group_by(year,AgeGroup) %>% summarise(ncase=n())

AG=10 # Number of age groups considered
max_A=49 # Max age considered in the case (note, does not have to be max age. The oldest cases may provide little info about transmission, as we
## are primarily interested in FOI over the period that serosurvey participants were alive)
T_case=8 # Time span of case data


# T_case=14 # Time span of case data






# Construct case data to input to model
case_data = matrix(nrow=AG,ncol=T_case)
agegroups = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
rownames(case_data) = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
colnames(case_data) = 2016:2023

# old cases
# colnames(case_data) = 2010:2023


for (r in 1:nrow(case_data)) {
      for (c in 1:ncol(case_data)) {
            case_data[r,c] = ifelse(length(caseage_by_year$ncase[caseage_by_year$year==colnames(case_data)[c] & caseage_by_year$AgeGroup==rownames(case_data)[r]])>0,
                                    caseage_by_year$ncase[caseage_by_year$year==colnames(case_data)[c] & caseage_by_year$AgeGroup==rownames(case_data)[r]],0)
      }
}
