

library(tidyverse)
library(tidybayes)


setwd("code/catalytic_model")


everyone <- read_rds(file="output/model_fits/20250603_2134_everyone.rds")
lifelong <- read_rds(file="output/model_fits/20250603_2138_lifelong.rds")
former <- read_rds(file="output/model_fits/20250603_2142_former.rds")
latter <- read_rds(file="output/model_fits/20250603_2146_latter.rds")


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



all_draws %>%
      mutate(year=2023-max(year)+year) %>%
      group_by(year,population)%>%
      median_qi() %>%
      filter(year>=2008)%>%
      mutate(Population=case_when(
            population=="everyone" ~ "All participants,\nno residency data",
            population=="lifelong" ~ "Only lifelong residents,\nno residency data",
            population=="latter" ~ "All participants,\nlatter location used",
            population=="former" ~ "All participants,\nformer location used"
      ))%>%
      ggplot(aes(x=factor(year),y=lambda*4,col=Population))+
      geom_point(aes(col=Population),
                 position = position_dodge(0.5)
                 )+
      geom_linerange(aes(col=Population,
                         ymin=.lower*4,ymax = .upper*4),
                     position = position_dodge(0.5)
                     )+
      xlab("Year")+
      ylab("Force of Infection")
      # geom_ribbon(aes(ymin=.lower*4,ymax = .upper*4),alpha=0.1,col=NA)+
      # facet_wrap(.~population,scales="free")
      

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


combo_fit <- read_rds(file="output/model_fits/20250603_2241_former_combo1718.rds")
combo_fit%>% 
      spread_draws(lambda[year]) %>%
      ungroup()%>%
      mutate(year=2023-max(year)+year) %>%
      mutate(newyear=case_when(
            year>2018 ~ as.character(year),
            year==2018 ~"2017-2018",
            year<2018 ~ as.character(year-1)
      )) %>%
      filter(year>=2009)%>%
      group_by(newyear)%>%
      median_qi() %>%
      ggplot(aes(x=newyear,y=lambda*4))+
      geom_point()+
      geom_linerange(aes(ymin=lambda.lower*4,ymax = lambda.upper*4))+
      scale_y_continuous(limits=c(0,1.5))+
      xlab("Year")+ylab("Force of Infection")+
      theme(axis.text.x = element_text(angle = 45,hjust=1))



lambdas <- combo_fit %>% spread_draws(lambda[t]) 

ages <- 1:49
times <-1:64

cumlambdas <- data.frame()

for(i in times){
      for(a in ages){
            
            print(i)
            print(a)
            
            tmp_cum_lambdas <- lambdas %>% 
                  filter(t<=i) %>% ### only include the years that have occurred before or during time i
                  filter(t>=(i-a)) %>% ### only include those years where the individuals are alive
                  group_by(.draw)%>%
                  summarise(cumlambda=sum(lambda)) %>%
                  mutate(age=a,
                         year=i
                  )
            
            
            cumlambdas <- bind_rows(cumlambdas,
                                    tmp_cum_lambdas
            )
            
            
            
            
            
      }
}




pop_groups<- cumlambdas %>%
      mutate(susc=exp(-4 * cumlambda),
             mono= 4*exp(-3*cumlambda)*(1-exp(-cumlambda)),
      )%>%
      mutate(the_rest=1-susc-mono)%>%
      group_by(.draw,year)%>%
      summarise(median=mean(cumlambda),
                susc=mean(susc),
                mono=mean(mono),
                the_rest=mean(the_rest)
      ) 


pop_groups %>%
      select(-median) %>%
      gather(type,value,-c(year,.draw)) %>%
      group_by(type,year)%>%
      summarise(median=median(value),
                lb=quantile(value,0.025),
                ub=quantile(value,0.975)
      ) %>%
      mutate(Type=factor(type,levels=c("susc", "mono", "the_rest"),labels=c("DENV-naïve","Primary", "Multiple")))%>%
      mutate(year=2023-max(year)+year) %>%
      mutate(newyear=case_when(
            year>2018 ~ as.character(year),
            year==2018 ~"2017-2018",
            year<2018 ~ as.character(year-1)
      )) %>%
      filter(year>=2009)%>%
      ggplot(aes(x=newyear,col=Type,fill=Type,group=Type))+
      geom_line(aes(y=median))+
      geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5,col=NA)+
      labs(
            #title = "Trends in Population Immunity Over Time",
            x = "Year",
            y = "Proportion of population by DENV infection History",
            color = "Infection History",
            fill = "Infection History"
      ) +
      theme_minimal()+
      theme(axis.text.x = element_text(angle = 45,hjust=1))



new_fig3<- cowplot::plot_grid(
      former$p1 +cowplot::theme_cowplot(),
      cowplot::plot_grid(former$p2 +cowplot::theme_cowplot(),
                         former$p3 +cowplot::theme_cowplot(),nrow = 1
      ), ncol=1)


new_fig3
# 
# 
# #### compare model fits
# everyone <- read_rds(file="output/20250527_2059_everyone.rds")
# lifelong <- read_rds(file="output/20250527_212_lifelong.rds")
# former <- read_rds(file="output/20250527_216_former.rds")
# latter <- read_rds(file="output/20250527_2111_latter.rds")
# 
# 
# everyone_prev <- read_rds(file="output/20250519_1624_everyone.rds")
# lifelong_prev <- read_rds(file="output/20250519_1628_lifelong.rds")
# former_prev <- read_rds(file="output/20250519_1633_former.rds")
# latter_prev <- read_rds(file="output/20250519_1639_latter.rds")
# 
# 
# 
# 
# everyone$fit %>% shinystan::launch_shinystan()
# 
# former$combined_plot
# 
# 
# 
# all_draws <- bind_rows(bind_rows(
#       former$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="former"),
#       latter$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="latter"),
#       
#       everyone$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="everyone"),
#       lifelong$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="lifelong")
#       
#             )%>% mutate(model="new"),
#       bind_rows(
#       former_prev$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="former"),
#       latter_prev$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="latter"),
#       
#       everyone_prev$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="everyone"),
#       lifelong_prev$fit %>% 
#             spread_draws(lambda[year]) %>%
#             mutate(population="lifelong")
#       
#       )%>% mutate(model="previous")
#       ) %>%
#       ungroup()
# 
# 
# 
# 
# 
# 
# 
# all_draws %>%
#       # filter(population=="former")%>%
#       mutate(year=2023-max(year)+year) %>%
#       # mutate(lambda4=lambda*4)%>%
#       group_by(year,population,model)%>%
#       median_qi() %>%
#       filter(year>=2005)%>%
#       ggplot(aes(x=year,y=lambda*4,col=model,fill=model))+
#       geom_line()+
#       geom_ribbon(aes(ymin=.lower*4,ymax = .upper*4),alpha=0.1,col=NA)+
#       facet_wrap(.~population,scales="free")
# 
# 

