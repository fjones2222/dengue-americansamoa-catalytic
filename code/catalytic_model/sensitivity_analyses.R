

library(tidyverse)
library(tidybayes)


setwd("code/catalytic_model")


everyone <- read_rds(file="output/20250528_1832_everyone.rds")
lifelong <- read_rds(file="output/20250528_1836_lifelong.rds")
former <- read_rds(file="output/20250528_1840_former.rds")
latter <- read_rds(file="output/20250528_1844_latter.rds")


everyone$fit %>% rhat() %>% hist()
lifelong$fit %>% rhat() %>% hist()
former$fit %>% rhat() %>% hist()
latter$fit %>% rhat() %>% hist()

everyone$fit %>% shinystan::launch_shinystan()

former$combined_plot



all_draws <- bind_rows(
      former$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="former"),
      latter$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="latter"),
      
      everyone$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="everyone"),
      lifelong$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="lifelong"),
      
) %>%
      ungroup()


all_draws %>%
      # filter(population=="former")%>%
      mutate(year=2023-max(year)+year) %>%
      # mutate(lambda4=lambda*4)%>%
      group_by(year,population)%>%
      median_qi() %>%
      filter(year>=2005)%>%
      ggplot(aes(x=year,y=lambda*4,col=population,fill=population))+
      geom_line()+
      geom_ribbon(aes(ymin=.lower*4,ymax = .upper*4),alpha=0.1,col=NA)+
      facet_wrap(.~population,scales="free")
      

former$fit %>%
      spread_draws(cum_lambda[i,j]) %>% summary()

former$p1

former$fit %>% spread_draws(lambda_temp[i],lambda[i]) %>% 
      mutate(year=2023-64+i)%>%
      filter(year %in% 2006:2016)%>%
      ggplot(aes(x=lambda,y=lambda_temp))+
      geom_point()+
      facet_wrap(.~year)

former_4sero <- read_rds(file="output/20250507_1123_former_4serotypes.rds")
former_3sero <- read_rds(file="output/20250507_1117_former_3serotypes.rds")
former_2sero <- read_rds(file="output/20250507_1110_former_2serotypes.rds")



sero_draws <- bind_rows(
      former_4sero$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(serotypes=4),
      former_3sero$fit %>% 
            spread_draws(lambda[year])  %>%
            mutate(serotypes=3),
      former_2sero$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(serotypes=2)
) %>%
      ungroup()




sero_draws %>%
      mutate(year=2023-max(year)+year) %>%
      # mutate(lambda4=lambda*4)%>%
      group_by(year,serotypes)%>%
      median_qi() %>%
      filter(year>=2010)%>%
      ggplot(aes(x=year,y=lambda*serotypes,col=factor(serotypes),fill=factor(serotypes)))+
      geom_line()+
      geom_ribbon(aes(ymin=.lower*serotypes,ymax = .upper*serotypes),alpha=0.1,col=NA)



read_rds(file="output/20250507_226_former_combo1718.rds")%>% 
      spread_draws(lambda[year]) %>%
      ungroup()%>%
      mutate(year=2023-max(year)+year) %>%
      group_by(year)%>%
      median_qi() %>%
      filter(year>=2011)%>%
      ggplot(aes(x=year,y=lambda*4))+
      geom_line()+
      geom_ribbon(aes(ymin=.lower*4,ymax = .upper*4),alpha=0.1,col=NA)+
      scale_y_continuous(limits=c(0,1.5))








#### compare model fits
everyone <- read_rds(file="output/20250527_2059_everyone.rds")
lifelong <- read_rds(file="output/20250527_212_lifelong.rds")
former <- read_rds(file="output/20250527_216_former.rds")
latter <- read_rds(file="output/20250527_2111_latter.rds")


everyone_prev <- read_rds(file="output/20250519_1624_everyone.rds")
lifelong_prev <- read_rds(file="output/20250519_1628_lifelong.rds")
former_prev <- read_rds(file="output/20250519_1633_former.rds")
latter_prev <- read_rds(file="output/20250519_1639_latter.rds")




everyone$fit %>% shinystan::launch_shinystan()

former$combined_plot



all_draws <- bind_rows(bind_rows(
      former$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="former"),
      latter$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="latter"),
      
      everyone$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="everyone"),
      lifelong$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="lifelong")
      
            )%>% mutate(model="new"),
      bind_rows(
      former_prev$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="former"),
      latter_prev$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="latter"),
      
      everyone_prev$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="everyone"),
      lifelong_prev$fit %>% 
            spread_draws(lambda[year]) %>%
            mutate(population="lifelong")
      
      )%>% mutate(model="previous")
      ) %>%
      ungroup()







all_draws %>%
      # filter(population=="former")%>%
      mutate(year=2023-max(year)+year) %>%
      # mutate(lambda4=lambda*4)%>%
      group_by(year,population,model)%>%
      median_qi() %>%
      filter(year>=2005)%>%
      ggplot(aes(x=year,y=lambda*4,col=model,fill=model))+
      geom_line()+
      geom_ribbon(aes(ymin=.lower*4,ymax = .upper*4),alpha=0.1,col=NA)+
      facet_wrap(.~population,scales="free")



