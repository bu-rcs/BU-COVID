## ---------------------------
##
## Generate classroom networks
##
## Authors: Wenrui Li, Katia Bulekova 
##          Boston University
##
## Date Created: 2020-07-31
##
## Email: wenruili@bu.edu
##
## ---------------------------

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# covasim output csv 
dir_csv <- "../Data/results/no_interventions/sim_results_no_interventions.csv"

# population size
n_pop <- 3750

# get cumulative infections, cumulative recoveries, cumulative deaths
sim_results <- read_csv(dir_csv) %>% 
  select(sim_num,dates, cum_infections, cum_recoveries, cum_deaths) %>%
  rename(Infections="cum_infections", Recoveries="cum_recoveries", Deaths="cum_deaths") %>%
  gather(Type, counts, Infections:Deaths, factor_key = TRUE ) %>%
  group_by(dates,Type) %>% 
  summarise(means=mean(counts),
            Low=quantile(counts,probs = 0.1),
            High=quantile(counts,probs = 0.9))  

# make plot
ggplot(sim_results, aes(as.Date(dates), y = means/n_pop,color=Type) ) +
  geom_ribbon( aes(ymin = Low/n_pop, 
                   ymax = High/n_pop, 
                   fill = Type, 
                   color = NULL), 
               alpha = .15) +
  geom_line( aes(y = means/n_pop), size = 1) + 
  xlab("") +
  ylab("Percent") +
  ggtitle("Cumulative percentages of model population") +
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d") +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),legend.position = c(0.1, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())


