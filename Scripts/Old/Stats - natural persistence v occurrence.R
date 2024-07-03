
#-----------------------------------------------------------------------------.
# DESCRIPTION: Figure 1c persistence & occurrence proportions in natural conditions with neighbors
#-----------------------------------------------------------------------------.

###############################
# Load packages and functions
###############################

# library(readr)
# library(tidyverse)
# library(truncnorm)
# library(stringr)
# library(BiodiversityR)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(visreg)
# library(Rmisc)
# library(patchwork)
# library(rgeos)

source("Scripts/Source - MAIN fitnessdata.R")

ci_prop <- function(level = 0.975, n, p) qt(level,df=n-1)*sqrt(p*(1-p)/n)

{
  ########################
  # DATA WRANGLING
  ########################
  
  # calculate proportion
  plotlev$persistence <- as.factor(plotlev$persistence)
  plotlev$occurrence <- as.factor(plotlev$occurrence)
  plotlev$names <- as.factor(plotlev$names)
  
  #occurrence
  fig1A_prop1 <- plotlev %>%
    filter(treatment %in% 'B') %>%
    select(names, occurrence, species) %>%
    distinct() %>%
    dplyr::group_by(names) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(names, occurrence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  N_occ <- length(unique(fig1A_prop1$names)) # use for CI calc
  fig1A_prop1_yes <-fig1A_prop1 %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("occurrence")) %>%
    select(-occurrence, -replicates)
  fig1A_prop1_no <- fig1A_prop1 %>% 
    filter(occurrence %in% "no") %>% # idea is to keep zeros, so find where occurrence no = 1, and make the prop zero, then rbind to persistence yes dataframe
    mutate(treatment = c("occurrence")) %>%
    filter(prop == 1) %>%
    select(-occurrence, -replicates)
  fig1A_prop1_no$prop[fig1A_prop1_no$prop == 1] <- 0    
  fig1A_prop1 <- rbind(fig1A_prop1_yes,fig1A_prop1_no)  
  
  #persitence
  fig1A_prop2 <- plotlev %>%
    filter(treatment %in% 'B') %>%
    select(names, persistence, species) %>%
    distinct() %>%
    dplyr::group_by(names) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(names, persistence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  N_per <- length(unique(fig1A_prop2$names)) # use for CI calc
  fig1A_prop2_yes <-fig1A_prop2 %>% # now filter for 'yes' 
    filter(persistence %in% 'yes') %>%
    mutate(treatment = c("persistence")) %>%
    select(-persistence, -replicates)
  fig1A_prop2_no <- fig1A_prop2 %>% 
    filter(persistence %in% "no") %>% # idea is to keep zeros, so find where persistence no = 1, and make the prop zero, then rbind to peresistence yes dataframe
    filter(prop == 1) %>%
    mutate(treatment = c("persistence")) %>%
    select(-persistence, -replicates)
  fig1A_prop2_no$prop[fig1A_prop2_no$prop == 1] <- 0    
  fig1A_prop2 <- rbind(fig1A_prop2_yes,fig1A_prop2_no)  
  
  
  fig1A_prop <- rbind(fig1A_prop1,fig1A_prop2)
  
  fig1A_prop$colortreat <- NA
  fig1A_prop$colortreat[fig1A_prop$treatment == 'occurrence']<- "midnightblue" # was lightgrey for dark mode
  fig1A_prop$colortreat[fig1A_prop$treatment == 'persistence']<- "cornflowerblue"

  temp_key <- nkey %>%
    select(names, site, grid) %>%
    distinct()
  fig1A_prop <- left_join(fig1A_prop, temp_key, by = 'names') # join to account for random effects
  
  # correct variable types
  fig1A_prop$names<- as.factor(fig1A_prop$names)
  fig1A_prop$treatment<- as.factor(fig1A_prop$treatment)
  fig1A_prop$colortreat<- as.factor(fig1A_prop$colortreat)
  fig1A_prop$site<- as.factor(fig1A_prop$site)
  fig1A_prop$grid<- as.factor(fig1A_prop$grid)
  
  } 

{
########################
# Fit model
########################

  hist(fig1A_prop$prop)
  
  # fit abiotic and biotic separately as they should not add to 1 (independent treatments)
  mm1o <- glmmTMB::glmmTMB(prop ~ 1 + (1|grid/site), data = filter(fig1A_prop, treatment %in% 'occurrence')) 
  summary(mm1o)
  
  mm1p <-  glmmTMB::glmmTMB(prop ~ 1 + (1|grid/site), data = filter(fig1A_prop, treatment %in% 'persistence')) 
  summary(mm1p)
  
  z_o <- summary(mm1o)$coefficient$cond[1,1] # occurrence estimate
  z_p <- summary(mm1p)$coefficient$cond[1,1] # persistence estimate
  
  Anova(mm1o, treatment = 3)
  performance::r2(mm1o)
  
  Anova(mm1p, treatment = 3)
  performance::r2(mm1p)
  
  
  # check if actual mean is far off from model estimate
  fig1A_prop %>%
    dplyr::group_by(treatment) %>%
    dplyr::mutate(mean_prop = mean(prop)) %>%
    select(mean_prop, treatment) %>%
    distinct() # yes the same !!! Good 
  
  
  ########################
  # CONFIDENCE INTERVALS
  ########################
  
  prop <- c(z_o,z_p)
  low <- c(as.numeric(prop[1] - ci_prop(n=N_occ[1], p=prop[1])), as.numeric(prop[2]  - ci_prop(n=N_per[1], p=prop[2]))) 
  high <- c(as.numeric(prop[1] + ci_prop(n=N_occ[1], p=prop[1])), as.numeric(prop[2]  + ci_prop(n=N_per[1], p=prop[2]))) 
  treatment <- c("occurrence","suitability")
  summary_1A <- data.frame(treatment=treatment,prop = prop, low=low, high=high)
  summary_1A$treatment <- as.factor(summary_1A$treatment)
  
  # persistence is same as persistence B in next figure. good check.
}
