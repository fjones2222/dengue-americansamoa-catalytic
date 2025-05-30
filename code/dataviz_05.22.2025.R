#installpackages

pacman::p_load(tidyverse)
install.packages("readxl")
install.packages("cowplot")


#load library packages 
library(readxl)
library(ggplot2)
library(lubridate)
library(tidyverse)
library(tidybayes)


## library for plots###
library(ggplot2)
library(grid)  # For unit function
library(cowplot)

#load dataset from American Samoa Serosurvey

analytical_dataset_brief<-analytical_dataset %>%
      select (record_id, age, yob, final_test_result)


#load case data
Arbo <- read_excel("ArboNET.xlsx")
Arbo$age <- as.numeric(Arbo$age)
Arbo <- Arbo %>%
      mutate(
            case = 1,
            age_cat = ifelse(is.na(age), NA, ifelse(age < 18, "<18", ">=18"))
      )

summarized_da

#adding 0s from 2019-2023

####----Analysis ArboNET----####




#2010 serosurvey data
# Creating the data frame again
sero_2010 <- data.frame(
      Year = rep(2010, 4),
      Age_Group = c("18-25", "26-40", "41-53", "54-87"),
      Positive = c(179, 216, 182, 182),
      Total = c(201, 217, 187, 189)
)

#calculate the seroprevalence
sero_2010$Seroprevalence <- (sero_2010$Positive / sero_2010$Total) * 100

sero_2010 <- sero_2010 %>%
      select(Year, Age_Group, Positive, Seroprevalence)


#combine the datasets
combined_datasets_AS <- bind_rows(Arbo, analytical_dataset_brief, sero_2010)

# Data preparation
summarized_data <- combined_datasets_AS %>%
      group_by(year, age_cat) %>%
      summarize(case = n(), .groups = 'drop') %>%
      filter(!is.na(case) & !is.na(year) & !is.na(age_cat))

# Ensure years 2019–2021 are included with 0 cases if missing
fill_years <- data.frame(
      year = rep(2019:2021, each = 2),
      age_cat = rep(c("<18", ">=18"), times = 3),
      case = 0
)

# Add to summarized_data
summarized_data <- bind_rows(summarized_data, fill_years)


# Ensure max_case is calculated
max_case <- max(summarized_data$case, na.rm = TRUE)

# Aggregate yearly totals
yearly_totals <- aggregate(case ~ year, data = summarized_data, sum)

# Create base plot
p <- ggplot(summarized_data, aes(x = year, y = case)) +
      geom_bar(aes(fill = age_cat), stat = "identity", position = "stack") +
      geom_text(data = yearly_totals, 
                aes(x = year, y = case, label = case),
                vjust = -0.5,
                color = "black", 
                size = 3.5) +
      scale_fill_manual(values = c("<18" = "#de2d26", ">=18" = "#3182bd")) +
      labs(x = "Year", 
           y = "Number of Cases",
           caption = "Note: Seroprevalence data from national serosurveys conducted in 2010 and 2023") +
      theme_minimal() +
      theme(
            legend.position = c(0.8, 1),  
            legend.justification = c(0.5, 1),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            plot.margin = margin(20, 5, 5, 5),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(angle = 0, hjust = 0.5)
      ) +
      guides(fill = guide_legend(title = "Age Category"))

# Final annotated plot
plot1 <- p +
      scale_x_continuous(
            limits = c(2007, 2025),  # Extended to give 2023 room
            breaks = seq(2008, 2023, by = 1),
            expand = expansion(mult = c(0.01, 0.02))
      ) +
      
      # Serosurvey annotation boxes
      annotate("label", x = 2010.5, y = max_case * 0.50, 
               label = "2010 Serosurvey\n18–25: 89.0%\n26–40: 99.5%\nn=395",
               size = 4, color = "black", fill = "white", label.size = 0.5, 
               lineheight = 1.3, hjust = 0.5, vjust = 0.5) +
      annotate("label", x = 2023, y = max_case * 0.6, 
               label = "2023 Serosurvey\n7–16: 59%\nn=887",
               size = 4, color = "black", fill = "white", label.size = 0.5, 
               lineheight = 1.2, hjust = 0.5, vjust = 0.5) +
      
      # Dotted lines for outbreaks
      annotate("segment", 
               x = 2009, xend = 2009,
               y = 0, yend = max_case * 0.87,
               linetype = "dotted", color = "darkred", linewidth = 1) +
      annotate("segment", 
               x = 2015, xend = 2015,
               y = 0, yend = max_case * 0.96,
               linetype = "dotted", color = "darkred", linewidth = 1) +
      
      # Outbreak labels with two lines
      annotate("label", x = 2009, y = max_case * 0.89,
               label = "DENV-4\noutbreak?",
               size = 4.8, color = "black", fill = "white", label.size = 0.5,
               hjust = 0.5, lineheight = 1.1) +
      annotate("label", x = 2015, y = max_case * 0.97,
               label = "DENV-2\noutbreak?",
               size = 4.8, color = "black", fill = "white", label.size = 0.5,
               hjust = 0.5, lineheight = 1.1) +
      
      # Arrows from serosurvey boxes
      annotate("segment", 
               x = 2010, xend = 2010,
               y = max_case * 0.3, yend = 0,  
               arrow = arrow(length = unit(0.7, "cm")), 
               color = "black", linewidth = 1.2) +
      annotate("segment", 
               x = 2023, xend = 2023,
               y = max_case * 0.4, yend = 0, 
               arrow = arrow(length = unit(0.7, "cm")), 
               color = "black", linewidth = 1.2)

# Display the plot
plot1

ggsave("Figures/epicurve1.png", plot = plot1, width = 10, height = 6, bg = "white")


#Dataset 2: Residency Figures 
### with weighted data####
library(readxl)
wgt_sum <- read_xlsx("wgt_sum.xlsx")
wgt_sum <- wgt_sum %>%
      mutate(
            point_prop = point / 100,
            lower_prop = lower / 100,
            upper_prop = upper / 100
      )

f2a <- ggplot(wgt_sum, aes(x = birth_group, y = point_prop)) +
      geom_point(size = 3, color = "blue") +  # Point estimates
      geom_errorbar(aes(ymin = lower_prop, ymax = upper_prop), width = 0.2, color = "black") +  # Error bars for CIs
      scale_y_continuous(limits = c(0, 1), 
                         breaks = seq(0, 1, by = 0.2), 
                         labels = scales::percent) +
      labs(x = "Birthplace location group",
           y = "Prevalence of previous dengue infection") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

f2a

###weightedestimates###
f2a <- survey_data1%>%
      mutate(location_group = case_when(
            birthplace_clean == "Missing" ~ NA_character_,  # Use NA_character_ for character columns
            birthplace_clean == "American Samoa" ~ "American Samoa",
            birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            birthplace_clean %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            TRUE ~ birthplace_clean # Retain any other locations unchanged
      ))%>%
      group_by(location_group)%>%
      summarise(seropos=mean(`Final Test Result`== "P")) %>%
      ggplot(aes(x=location_group,y=seropos))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))



###weightedestimates###
f2a <- survey_data1%>%
      mutate(location_group = case_when(
            birthplace_clean == "Missing" ~ NA_character_,  # Use NA_character_ for character columns
            birthplace_clean == "American Samoa" ~ "American Samoa",
            birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            birthplace_clean %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            TRUE ~ birthplace_clean # Retain any other locations unchanged
      ))%>%
      group_by(location_group)%>%
      summarise(seropos=mean(`Final Test Result`== "P")) %>%
      ggplot(aes(x=location_group,y=seropos))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))

##create dot plots-include weighted estimates###
f2b <- analytical_long%>%
      mutate(birth_group = case_when(
            birthplace_clean == "Missing" ~ NA_character_,  # Use NA_character_ for character columns
            birthplace_clean == "American Samoa" ~ "American Samoa",
            birthplace_clean %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            birthplace_clean %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            TRUE ~ birthplace_clean # Retain any other locations unchanged
      ))%>%
      mutate(location_group = case_when(
            former_location == "Missing" ~ NA_character_,  # Use NA_character_ for character columns
            former_location == "American Samoa" ~ "American Samoa",
            former_location %in% c("Western Samoa", "Fiji", "Tuvalu", "Philippines", "Hawaii") ~ "Pacific Islands",
            former_location %in% c("New Zealand", "USA", "Australia") ~ "Non-endemic",
            TRUE ~ former_location # Retain any other locations unchanged
      ))%>%
      filter(birth_group=="Non-endemic")%>%
      left_join(select(analytical_dataset,
                       record_id=record_id,
                       final_test_result))%>%
      mutate(final_test_result=ifelse(final_test_result=="P",
                                        "Positive","Negative"
      ))%>%
      filter(!is.na(location_group))%>%
      ggplot(aes(x=factor(year), y=reorder(factor(record_id),yob)))+
      geom_tile(aes(fill=location_group),col="black")+
      geom_vline(
            xintercept = which(levels(factor(f2b$year)) == "2017"),
            color = "black",
            linewidth = 3  
      )+
      facet_grid(final_test_result~., scales="free_y", space = "free_y")+
      ylab("Serosurvey Participant ID")+
      xlab("Year")+
      scale_fill_brewer("Location", palette="Set2")+
      theme_cowplot()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.text.x = element_text(angle=45,hjust = 1)
      )
f2b


#cowplot
# combined plots
# Adjust f2a theme
f2a <- f2a + theme_bw() + 
      theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

# Adjust f2b theme and remove y-axis participant IDs
f2b <- f2b + theme_bw() + 
      theme(
            plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # Rotate x-axis labels
            axis.text.y = element_blank(),  # Remove participant IDs
            axis.ticks.y = element_blank()   # Remove y-axis ticks
      ) +
      ylab(NULL) 

# Combine the plots using cowplot
combined_plot <- cowplot::plot_grid(
      f2a, f2b,
      labels = c("A", "B"),
      align = 'h',  # Align horizontally
      axis = 'tb',  # Align both top and bottom
      nrow = 1,     # Single row layout
      rel_widths = c(1, 1)  # Equal widths
)


library(tidybayes)
print(combined_plot)

#Research Question 1: 1.Is there evidence of local transmission of dengue in American Samoa during years without cases reported?
modelobject<- read_rds("code/catalytic_model/20250519_1633_former.rds")
##modelfit<-read_rds("code/catalytic_model/as_denguefoi_fixeff_casess_Tlambda_17_mi_3_covs__ns_4_case2010to2022.rds")
modelfit<-modelobject$fit
draws<-modelfit%>%
      spread_draws(lambda[year])
summary(draws)
#mutate quantiles
#group by year
#year=2023+63
filter(draws, year>=44)%>%
      ggplot(aes(x=year, y=lambda*4, group=year))+
      geom_boxplot()

#refining this based on the reccs 
draws1 <- modelfit %>%
      spread_draws(lambda[year]) %>%
      group_by(year) %>%
      summarize(
            median = median(lambda * 4),
            lower = quantile(lambda * 4, 0.25),
            upper = quantile(lambda * 4, 0.75),
            .groups = 'drop'
      ) %>%
      mutate(actual_year = 2023 + (year - 64))%>%  # Convert index to actual years
      filter(actual_year >= 2008) 
      
 ###foi is draws 1     
      draws1
f4a<- ggplot(draws1, aes(x = actual_year)) +
      geom_point(aes(y = median)) +
      geom_linerange(aes(ymin = lower, ymax = upper)) +
      labs(x = "Year", y = "Force of infection") +
      scale_x_continuous(breaks = 2008:2023) +  # Show all years
      theme_minimal() +  # Using minimal theme
      theme(
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.background = element_blank(),   # Remove panel background
            plot.background = element_blank(),    # Remove plot background
            axis.line = element_line(color = "black")  # Add axis lines
      )

f4a

draws<-modelfit%>%spread_draws(alpha, lambda[t])

##### Research Question 2 Updated

former <- read_rds(file="code/catalytic_model/20250519_1633_former.rds")

lambdas <- former$fit %>% spread_draws(lambda[t]) 

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


f4b <-pop_groups %>%
      select(-median) %>%
      gather(type,value,-c(year,.draw)) %>%
      group_by(type,year)%>%
      summarise(median=median(value),
                lb=quantile(value,0.025),
                ub=quantile(value,0.975)
      ) %>%
      mutate(year=2023-64+year)%>%
      mutate(Type=factor(type,levels=c("susc", "mono", "the_rest"),labels=c("DENV-naïve","Primary", "Multiple")))%>%
      filter(year>2005)%>%
      ggplot(aes(x=year,col=type,fill=type))+
      geom_line(aes(y=median))+
      geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.5,col=NA)+
      labs(
            #title = "Trends in Population Immunity Over Time",
            x = "Year",
            y = "Proportion of population by DENV infection History",
            color = "Infection History",
            fill = "Infection History"
      ) +
      theme_minimal()



####Research Question 2: immunity before the previous outbreak and what is it now?

census <- read_csv("code/catalytic_model/Data/updated_2020census.csv") 
census_long <- census %>%
      mutate(ag1=as.numeric(substr(age_group,0,2)),
             ag2=ag1+1,
             ag3=ag2+1,
             ag4=ag3+1,
             ag5=ag4+1
      ) %>%
      mutate(pop=number/5) %>%
      select(ag1:pop) %>%
      gather(type,age,-pop) %>%
      select(age,pop) %>%
      arrange(age)
draws4 <- modelfit%>%
      spread_draws(susc[age, year], mono[age,year])%>%
      mutate(protct=1-mono-susc) %>%
      left_join(census_long)
t4b<-
      draws4 %>%
      group_by(year, .draw) %>%
      summarize(mono = mean(mono),
                protct = mean(protct),
                susc = mean(susc)) %>%
      gather(type, value, -c(year, .draw)) %>%
      group_by(year, type) %>%
      summarise(median = median(value),
                lb = quantile(value, 0.025),
                ub = quantile(value, 0.975))%>%
      mutate(Type=factor(type,levels=c("susc", "mono", "protct"),labels=c("DENV-naïve","Primary", "Multiple")))

f4b <- t4b%>%
      ggplot(aes(x = year, y = median, col = Type, fill = Type)) +
      geom_line() +
      geom_ribbon(alpha = 0.5, aes(ymin = lb, ymax = ub)) +
      scale_color_manual(values = c("gold","#F8766D","#8BEDA7")) +
      scale_fill_manual(values =  c("gold","#F8766D","#8BEDA7")) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.7)) +  
      scale_x_continuous(breaks = seq(1, 8, 1), 
                         labels = 2016:2023) +
      labs(x = "Year",
           y = "Proportion of population\nby DENV Infection History")+
      theme_minimal() +
      theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.line = element_line(color = "black")
      )

f4b
# Line plot with ribbon

#this is for year 2023 
draws3_2023 <- modelfit %>%
      spread_draws(mono[age, year]) %>%
      filter(year == 7) %>%
      group_by(age) %>%
      summarise(
            median = median(mono),
            lb = quantile(mono, 0.025),
            ub = quantile(mono, 0.975),
            .groups = 'drop'
      )

f4c <- ggplot(draws3_2023, aes(x = age, y = median)) +
      geom_line() +
      geom_ribbon(alpha = 0.5, aes(ymin = lb, ymax = ub)) +
      labs(x = "Age (years)",
           y = "Proportion of people by age with exactly\n one previous DENV infection") +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.5)) + 
      scale_x_continuous(limits = c(0, 50), breaks = seq(0, 50, 5)) +  # Set limits and breaks
      theme_minimal() +
      theme(
            panel.grid.minor = element_blank(),  # Remove minor grid lines
            panel.grid.major = element_blank(),  # Remove major grid lines
            panel.background = element_blank(),   # Remove panel background
            plot.background = element_blank(),    # Remove plot background
            axis.line = element_line(color = "black")  # Add axis lines
      )
f4c

#all plots

f4a

f4b

f4c

# Modify f4a to show all years but angled for readability
f4a <- f4a + 
      scale_x_continuous(breaks = seq(2010, 2023, by = 1)) +  # Show every year
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))  # Angle the year labels and adjust their position

# Arrange horizontally
figure4 <- cowplot::plot_grid(f4a, f4b, f4c,
                   labels = c("A", "B", "C"),
                   nrow = 3,          
                   label_size = 12,   
                   align = 'h')
           # Add more space between plots

figure4

ggsave(filename = "Figures/figure4.png",
       plot = figure4,
       height = 6.42,
       width = 4.34,
       scale=1.5
       
       )





###`10 year olds`

t4b10<-
      draws4 %>%
      filter(age==10)%>%
      group_by(year, .draw) %>%
      summarize(mono = mean(mono),
                protct = mean(protct),
                susc = mean(susc)) %>%
      gather(type, value, -c(year, .draw)) %>%
      group_by(year, type) %>%
      summarise(median = median(value),
                lb = quantile(value, 0.025),
                ub = quantile(value, 0.975))%>%
      mutate(Type=factor(type,levels=c("susc", "mono", "protct"),labels=c("DENV-naïve","Primary", "Multiple")))

t4b10%>%
      ggplot(aes(x = year, y = median, col = Type, fill = Type)) +
      geom_line() +
      geom_ribbon(alpha = 0.5, aes(ymin = lb, ymax = ub)) +
      scale_color_manual(values = c("gold","#F8766D","#8BEDA7")) +
      scale_fill_manual(values =  c("gold","#F8766D","#8BEDA7")) +
      scale_y_continuous(labels = scales::percent_format(), limits = c(0, 0.7)) +  
      scale_x_continuous(breaks = seq(1, 8, 1), 
                         labels = 2016:2023) +
      labs(x = "Year",
           y = "Proportion of population\nby DENV Infection History")+
      theme_minimal() +
      theme(
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(),
            axis.line = element_line(color = "black")
      )

t4b10

write.csv(t4b10, "C:\\Users\\sfz5\\OneDrive - CDC\\Dengue\\American Samoa Epi Aid\\AS_Residency_Analysis\\Data\\t4b10.csv")



f4c