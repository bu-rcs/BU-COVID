## ---------------------------
##
## Plot results for single run
##
## Authors: Wenrui Li, Katia Bulekova 
##          Boston University
##
## Date Created: 2020-08-03
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
dir_csv <- "../Data/results/no_interventions/sim_results_single_run.csv"

# population size
n_pop <- 3750

# get cumulative infections, cumulative recoveries, cumulative deaths
sim_cum <- read_csv(dir_csv) %>% 
  select(dates, cum_infections, cum_recoveries, cum_deaths) %>%
  rename(Infections="cum_infections", Recoveries="cum_recoveries", Deaths="cum_deaths") %>%
  gather(Type, counts, Infections:Deaths, factor_key = TRUE ) 


# get daily infections, cumulative recoveries, cumulative deaths
sim_daily <- read_csv(dir_csv) %>%  
  select(dates, new_infections, new_recoveries, new_deaths) %>%
  rename(Infections="new_infections", Recoveries="new_recoveries", Deaths="new_deaths") %>%
  gather(Type, counts, Infections:Deaths, factor_key = TRUE ) 

##### All the various things that can be plotted: 
#'cum_infections', 'cum_infectious', 'cum_tests', 'cum_diagnoses', 'cum_recoveries', 'cum_symptomatic', 
# 'cum_severe', 'cum_critical', 'cum_deaths', 'cum_quarantined', 'new_infections', 'new_infectious', 
# 'new_tests', 'new_diagnoses', 'new_recoveries', 'new_symptomatic', 'new_severe', 'new_critical', 
# 'new_deaths', 'new_quarantined', 'n_susceptible', 'n_exposed', 'n_infectious', 'n_symptomatic', 
# 'n_severe', 'n_critical', 'n_diagnosed', 'n_quarantined', 'prevalence', 'incidence', 'r_eff', 
# 'doubling_time', 'test_yield', 't', 'date'

# Remark     
# The prefix "new" is used for flow variables, i.e. counting new events (infections/deaths/recoveries) on each timestep
# The prefix "n" is used for stock variables, i.e. counting the total number in any given state (sus/inf/rec/etc) on any particular timestep
# The prefix "cum" is used for cumulative variables, i.e. counting the total number that have ever been in a given state at some point in the sim

# make plot
ggplot(sim_cum, aes(as.Date(dates), y = counts/n_pop,color=Type) ) +
  geom_line( aes(y = counts/n_pop), size = 1) + 
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

ggplot(sim_daily, aes(as.Date(dates), y = counts/n_pop,color=Type) ) +
  geom_line( aes(y = counts/n_pop), size = 1) + 
  xlab("") +
  ylab("Percent") +
  ggtitle("Daily percentages of model population") +
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d") +
  scale_y_continuous(labels=scales::percent) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(),legend.position = c(0.1, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())