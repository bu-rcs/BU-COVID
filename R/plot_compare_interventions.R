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


# population size
n_pop <- 3750

# covasim output file for different intervention
hdir <- "../Data/results/"
name_csv <- c("no_interventions/sim_results_no_interventions.csv",
              "social_distancing/sim_results_social_distancing.csv",
              "testing/sim_results_testing.csv")

n_int <- length(name_csv) -1
sim_results <- vector("list",n_int + 2)
int_seq <- c(" No intervention",
             "+ Classroom masks/distancing",
             "+ Testing",
             "Purely exogenous (~1 per day) ")
int_seq <- factor(int_seq, levels = int_seq)

# get cumulative infections 
for (i in 1:(n_int + 1)) {
  dir_csv <- paste0(hdir,name_csv[i])
  sim_results[[i]] <-  read_csv(dir_csv) %>% 
    select(sim_num, dates, cum_infections) %>% 
    group_by(dates) %>% 
    summarise(means=mean(cum_infections), Low=quantile(cum_infections, probs = 0.1),
              High=quantile(cum_infections, probs = 0.9)) %>%
    mutate(Intervention = int_seq[i]) 
}

# Purely exogenous (~1 per day)
n_poi <- 1
n_sum <- 1000
n_day <- length(sim_results[[n_int+1]]$dates)
set.seed(1)
infect_seq100 <- matrix(rpois(n_day*n_sum, n_poi),nrow=n_sum)
infect_cum100 <- t(apply(infect_seq100, 1, cumsum))
sim_results[[n_int + 2]] <- data.frame(dates=sim_results[[n_int]]$dates,
                               means=colMeans(infect_cum100),Low=apply(infect_cum100,2,quantile,probs=0.1),
                               High=apply(infect_cum100,2,quantile,probs=0.9),Intervention = int_seq[n_int + 2])


data_long <- Reduce(function(x, y) merge(x, y, all=TRUE), sim_results)  

# make plot
ggplot(data_long, aes(as.Date(dates), y = means/n_pop,color=Intervention) ) +
  geom_ribbon( aes(ymin = Low/n_pop, 
                   ymax = High/n_pop, 
                   fill = Intervention, 
                   color = NULL), 
               alpha = .15) +
  geom_line( aes(y = means/n_pop), size = 1) + 
  xlab("") +
  ylab("") +
  ggtitle("Cumulative infection")+ ylab("Percent Infected") +
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d") +
  scale_y_continuous(labels=scales::percent,limits=c(0,1)) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank(), legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())
