## Population data comes from US census data
popdata = read_excel('Data/2020 census.xlsx')
popdata = popdata %>% 
      mutate(AgeGroup = gsub(' years','',AgeGroup)) %>% 
      filter(AgeGroup %in% agegroups) %>% mutate(npop=Number)

pop_data = matrix(nrow=AG,ncol=T_case)
rownames(pop_data) = c("0 to 4","5 to 9","10 to 14","15 to 19","20 to 24","25 to 29","30 to 34","35 to 39","40 to 44","45 to 49")
colnames(pop_data) = 2016:2023

# old cases
# colnames(case_data) = 2010:2023



for (r in 1:nrow(pop_data)) {
      for (c in 1:ncol(pop_data)) {
            pop_data[r,c] = popdata$npop[popdata$AgeGroup==rownames(pop_data)[r]]
      }
}