
#------------------------------------------------------------------------------
# DESCRIPTION: Script corresponds to figure 2 showing how persistence and occurrence
# proportions for both neighbor treatments across discrete scales 
#------------------------------------------------------------------------------

##############################
# Load packages and functions
##############################

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

#-------------------------------
# Block w/ & w/out neighbors
#-------------------------------

########################
# Data wrangling
########################
{
  # Plot level persistence
  
  fig2A_per <- plotlev %>%
    select(names, persistence, treatment, species) %>%
    distinct() %>%
    dplyr::group_by(names, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(names, treatment, persistence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2A_per_yes <- fig2A_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2A_per_no <- fig2A_per %>% 
    filter(persistence %in% "no") %>% # idea is to keep zeros, so find where persistence no = 1, and make the prop zero, then rbind to peresistence yes dataframe
    filter(prop == 1) %>%
    select(-persistence)
  fig2A_per_no$prop[fig2A_per_no$prop == 1] <- 0    
  fig2A_per <- rbind(fig2A_per_yes,fig2A_per_no)  
  
  temp_key <- nkey %>%
    select(names, site, grid) %>%
    distinct()
  fig2A_per <- left_join(fig2A_per, temp_key, by = 'names') # join to account for random effects
  
  # Plot level occurrence
  
  fig2A_occ <- plotlev %>%
    filter(treatment %in% 'B') %>%
    select(names, occurrence, species, grid, site) %>%
    distinct() %>%
    dplyr::group_by(names) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(names, occurrence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2A_occ_yes <-fig2A_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2A_occ_no <- fig2A_occ %>% 
    filter(occurrence %in% "no") %>% # idea is to keep zeros, so find where occurrence no = 1, and make the prop zero, then rbind to persistence yes dataframe
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2A_occ_no$prop[fig2A_occ_no$prop == 1] <- 0    
  fig2A_occ <- rbind(fig2A_occ_yes,fig2A_occ_no)  

  # put persistence and occurrence raw data together for figure and models
  fig_fig_2A <- rbind(fig2A_occ,fig2A_per)
  
  # correct variable types
  fig_fig_2A$names<- as.factor(fig_fig_2A$names)
  fig_fig_2A$treatment<- as.factor(fig_fig_2A$treatment)
  fig_fig_2A$site<- as.factor(fig_fig_2A$site)
  fig_fig_2A$grid<- as.factor(fig_fig_2A$grid)
  
  fig_fig_2A$treatment <- case_match(fig_fig_2A$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2A$treatment)
  
}

########################
# Fit models
########################

mplot <- glmmTMB::glmmTMB(prop ~ treatment + (1|names) + (1|site) + (1|site:grid), 
                                                  weights = replicates,
                                                  data = fig_fig_2A,
                                                  family = binomial())

############################
# Set up viz and table data
############################

# ESTIMATES

vis.p <- ggpredict(mplot, 
                   terms = c("treatment"), 
                   type = "fe", allow.new.levels=TRUE)

vis.p$scale <- c("plot")

# ANOVAS

a <-Anova(mplot, type = 3)

p.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c(rep('plot',time = 2)),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type")
                            )

# add pairwise comparison table

em <- emmeans(mplot, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# occurrence                    0.559 0.0272 Inf     0.505     0.611
# suitability with neighbors    0.186 0.0177 Inf     0.154     0.223
# suitability without neighbors 0.327 0.0250 Inf     0.280     0.377

mplot.contrast <- as.data.frame(pairs(em))
pairs(em)
# contrast                                                   odds.ratio     SE  df null z.ratio p.value
# occurrence / suitability with neighbors                         5.534 0.4534 Inf    1  20.881  <.0001
# occurrence / suitability without neighbors                      2.609 0.2002 Inf    1  12.492  <.0001
# suitability with neighbors / suitability without neighbors      0.471 0.0404 Inf    1  -8.776  <.0001

#-------------------------------
# Grid w/ & w/out neighbors
#-------------------------------

#################
# Data wrangling
#################

{
  # Grid level persistence
  
  fig2C_per <- gridlev %>%
    select(site, grid, persistence, treatment, species) %>%
    distinct() %>%
    dplyr::group_by(grid, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(grid, treatment, persistence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2C_per_yes <- fig2C_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2C_per_no <- fig2C_per %>% 
    filter(persistence %in% "no") %>% # idea is to keep zeros, so find where persistence no = 1, and make the prop zero, then rbind to peresistence yes dataframe
    filter(prop == 1) %>%
    select(-persistence)
  fig2C_per_no$prop[fig2C_per_no$prop == 1] <- 0    
  fig2C_per <- rbind(fig2C_per_yes,fig2C_per_no)  
  
  # Grid level occurrence
  
  fig2C_occ <- gridlev %>%
    filter(treatment %in% 'B') %>%
    select(occurrence, species, grid, site) %>%
    distinct() %>%
    dplyr::group_by(grid) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(grid, occurrence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2C_occ_yes <-fig2C_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2C_occ_no <- fig2C_occ %>% 
    filter(occurrence %in% "no") %>% # idea is to keep zeros, so find where occurrence no = 1, and make the prop zero, then rbind to persistence yes dataframe
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2C_occ_no$prop[fig2C_occ_no$prop == 1] <- 0    
  fig2C_occ <- rbind(fig2C_occ_yes,fig2C_occ_no)  
 
  # together for figure and model
  fig_fig_2C <- rbind(fig2C_occ,fig2C_per)
  
  # correct variable types
  fig_fig_2C$treatment<- as.factor(fig_fig_2C$treatment)
  fig_fig_2C$site<- as.factor(fig_fig_2C$site)
  fig_fig_2C$grid<- as.factor(fig_fig_2C$grid)
  
  fig_fig_2C$treatment <- case_match(fig_fig_2C$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2C$treatment)
  
}


########################
# Fit models
########################

gplot <- glmmTMB::glmmTMB(prop ~ treatment + (1|site) + (1|site:grid), 
                          weights = replicates,
                          data = fig_fig_2C,
                          family = binomial())

############################
# Set up viz and table data
############################

# ESTIMATES

vis.g <- ggpredict(gplot, 
                   terms = c("treatment"), 
                   type = "fe", allow.new.levels=TRUE)

vis.g$scale <- c("grid")

# ANOVAS

a <-Anova(gplot, type = 3)

g.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                      scale = c(rep('grid',time = 2)),
                      Df =c(a$Df[1],a$Df[2]),
                      P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                      predictor = c("Intercept","Data type")
)

# add pairwise comparison table

em <- emmeans(gplot, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# occurrence                    0.952 0.0279 Inf    0.8567     0.985
# suitability with neighbors    0.147 0.0568 Inf    0.0659     0.295
# suitability without neighbors 0.523 0.0974 Inf    0.3374     0.702

gplot.contrast <- as.data.frame(pairs(em))
pairs(em)
# contrast                                                   odds.ratio     SE  df null z.ratio p.value
# occurrence / suitability with neighbors                       116.650 77.4172 Inf    1   7.171  <.0001
# occurrence / suitability without neighbors                     18.309 11.0153 Inf    1   4.833  <.0001
# suitability with neighbors / suitability without neighbors      0.157  0.0667 Inf    1  -4.358  <.0001


#-------------------------------
# Site w/ & w/out neighbors
#-------------------------------

#################
# Data wrangling
#################

# Site level persistence
{
  fig2E_per <- sitelev %>%
    select(site,persistence, treatment, species) %>%
    distinct() %>%
    dplyr::group_by(site, treatment) %>% # find number of species replicates in each plot
    dplyr::mutate(replicates = n()) %>%
    group_by(site, treatment, persistence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    select(-species) %>%
    distinct() %>%
    ungroup()
  
  fig2E_per_yes <- fig2E_per %>% 
    filter(persistence %in% "yes")%>%
    select(-persistence) # filter for yes, as 1 - proportion persisting = proportion not persisting
  fig2E_per_no <- fig2E_per %>% 
    filter(persistence %in% "no") %>% # idea is to keep zeros, so find where persistence no = 1, and make the prop zero, then rbind to peresistence yes dataframe
    filter(prop == 1) %>%
    select(-persistence)
  fig2E_per_no$prop[fig2E_per_no$prop == 1] <- 0    
  fig2E_per <- rbind(fig2E_per_yes,fig2E_per_no)  
  
  # Site level occurrence
  
  fig2E_occ <- sitelev %>%
    filter(treatment %in% 'B') %>%
    select(site, occurrence, species) %>%
    distinct() %>%
    dplyr::group_by(site) %>%
    dplyr::mutate(replicates = n()) %>%
    dplyr::group_by(site, occurrence) %>%
    dplyr::mutate(prop = n()/replicates) %>%
    ungroup() %>%
    select(-species) %>%
    distinct()
  
  fig2E_occ_yes <-fig2E_occ %>% # now filter for 'yes' 
    filter(occurrence %in% 'yes') %>%
    mutate(treatment = c("O")) %>%
    select(-occurrence)
  fig2E_occ_no <- fig2E_occ %>% 
    filter(occurrence %in% "no") %>% # idea is to keep zeros, so find where occurrence no = 1, and make the prop zero, then rbind to persistence yes dataframe
    mutate(treatment = c("O")) %>%
    filter(prop == 1) %>%
    select(-occurrence)
  fig2E_occ_no$prop[fig2E_occ_no$prop == 1] <- 0    
  fig2E_occ <- rbind(fig2E_occ_yes,fig2E_occ_no)
 
  # pull together for plotting and model
  fig_fig_2E <- rbind(fig2E_occ,fig2E_per)
  
  # correct variable types
  fig_fig_2E$treatment<- as.factor(fig_fig_2E$treatment)
  fig_fig_2E$site<- as.factor(fig_fig_2E$site)
  
  fig_fig_2E$treatment <- case_match(fig_fig_2E$treatment, 
                                     "A" ~ "suitability without neighbors", 
                                     "B" ~ "suitability with neighbors", 
                                     "O" ~ "occupancy", 
                                     .default = fig_fig_2E$treatment)
}

########################
# Fit models
########################

splot <- glmmTMB::glmmTMB(prop ~ treatment + (1|site), 
                          weights = replicates,
                          data = fig_fig_2E,
                          family = binomial())

############################
# Set up viz and table data
############################

# ESTIMATES

vis.s <- ggpredict(splot, 
                   terms = c("treatment"), 
                   type = "fe", allow.new.levels=TRUE)

vis.s$scale <- c("site")

# ANOVAS

a <-Anova(splot, type = 3)

s.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                      scale = c(rep('site',time = 2)),
                      Df =c(a$Df[1],a$Df[2]),
                      P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                      predictor = c("Intercept","Data type")
)

# add pairwise comparison table

em <- emmeans(splot, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# occurrence                    1.0000000 0.00000184 Inf 0.0000000 1.0000000
# suitability with neighbors    0.1666667 0.07607257 Inf 0.0639867 0.3691369
# suitability without neighbors 0.5833333 0.10063459 Inf 0.3834264 0.7591408

splot.contrast <- as.data.frame(pairs(em))
pairs(em)
# contrast                                                   odds.ratio     SE  df null z.ratio p.value
# occurrence / suitability with neighbors                      6.18e+10 1.4e+15 Inf    1   0.001  1.0000
# occurrence / suitability without neighbors                   8.83e+09 2.0e+14 Inf    1   0.001  1.0000
# suitability with neighbors / suitability without neighbors   0.00e+00 0.0e+00 Inf    1  -2.834  0.0128

###############################################
# used for figure and tables

summary <- rbind(vis.p,vis.g,vis.s)
names(summary)
summary <- as.data.frame(summary)
###############################################


rm(fig2A_per_yes,fig2A_per_no,fig2A_occ_yes,fig2A_occ_no, 
   fig2C_per_yes,fig2C_per_no,fig2C_occ_yes,fig2C_occ_no, 
   fig2E_per_yes,fig2E_per_no,fig2E_occ_yes,fig2E_occ_no)

#-------- the end