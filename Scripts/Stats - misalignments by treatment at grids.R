
#-----------------------------------------------------------------------------.
# DESCRIPTION: Fitting a glms to grid level data (n = 18) for both neighbor 
# treatments, to find the esitmated proportions of each (mis)alignment. 
# Then, calculate confidence intervals around these proportions, and compare 
# proportions to assess significance using emmeans
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
# library(mclogit)
# library(emmeans)

#source("Scripts/Source - MAIN fitnessdata.R")

#--------------------------------------------------------------------------.
# META-POP SCALE (GRIDS)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2D_dat <- gridlev %>%
  select(grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(grid, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(grid, treatment, contingency) %>%
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
# all NAs = 0 because I standardized by replicate
fig2D_dat[is.na(fig2D_dat)] <- 0
fig2D_dat <- pivot_longer(fig2D_dat, cols = 5:8, names_to = "contingency", values_to = "prop")

fig2D_dat$contingency<-as.factor(fig2D_dat$contingency)
fig2D_dat$site<-as.factor(fig2D_dat$site)
fig2D_dat$grid<-as.factor(fig2D_dat$grid)
fig2D_dat$treatment<-as.factor(fig2D_dat$treatment)

########################
## Fit models 
########################

################
# sinks
m.me.g <- glmmTMB(prop ~ treatment + (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'ME'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.me.g)
testDispersion(m.me.g) 
plot(fitted(m.me.g), residuals(m.me.g))
hist(residuals(m.me.g)) 
loc_sim_ouput <- simulateResiduals(m.me.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.me.g)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.1202     0.3049  -0.394    0.694    
# treatmentB    1.5661     0.3967   3.948 7.88e-05 ***

a <- Anova(m.me.g, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)  0.1553  1     0.6935    
# treatment   15.5871  1  7.879e-05 ***

# save anova values for table:
grid.me.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('sink','sink'))

em <- emmeans(m.me.g, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.470 0.0760 Inf     0.328     0.617
# B         0.809 0.0554 Inf     0.677     0.896

grid.me.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
#  A / B         0.209 0.0828 Inf    1  -3.948  0.0001

################
# d. limitation
m.dl.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'DL'), 
                  weights = replicates, 
                  family = binomial(),
                  ziformula=~1, # adding a zero inflation term fixes convergence
                  control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS"))) 

# Check model fit 
testZeroInflation(m.dl.g)
testDispersion(m.dl.g) 
plot(fitted(m.dl.g), residuals(m.dl.g))
hist(residuals(m.dl.g)) 
loc_sim_ouput <- simulateResiduals(m.dl.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.dl.g)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -3.5264     0.7174  -4.915 8.87e-07 ***
# treatmentB   -0.7077     1.2366  -0.572    0.567

a <- Anova(m.dl.g, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 24.1599  1  8.866e-07 ***
# treatment    0.3276  1     0.5671   

# save anova values for table:
grid.dl.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('dispersal limitation','dispersal limitation'))


em <- emmeans(m.dl.g, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.0286 0.0199 Inf   0.00716    0.1071
# B         0.0143 0.0142 Inf   0.00201    0.0945

grid.dl.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B          2.03 2.51 Inf    1   0.572  0.5671

################
# a. present
m.ap.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'SS_y'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.ap.g)
testDispersion(m.ap.g) 
plot(fitted(m.ap.g), residuals(m.ap.g))
hist(residuals(m.ap.g)) 
loc_sim_ouput <- simulateResiduals(m.ap.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.ap.g)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -0.3490     0.3864  -0.903    0.366    
# treatmentB   -1.8976     0.4717  -4.023 5.76e-05 ***

a <- Anova(m.ap.g, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)  0.816  1     0.3664    
# treatment   16.181  1  5.758e-05 ***

# save anova values for table:
grid.ap.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned present','aligned present'))

em <- emmeans(m.ap.g, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.4136 0.0937 Inf    0.2486     0.601
# B         0.0956 0.0433 Inf    0.0381     0.220

grid.ap.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B          6.67 3.15 Inf    1   4.023  0.0001

################
# a. absent
m.aa.g <- glmmTMB(prop ~ treatment +  (1 | site) + (1| grid:site), 
                  data = filter(fig2D_dat, contingency %in% 'SS_n'), 
                  weights = replicates, 
                  family = binomial())

# Check model fit 
testZeroInflation(m.aa.g)
testDispersion(m.aa.g) 
plot(fitted(m.aa.g), residuals(m.aa.g))
hist(residuals(m.aa.g)) 
loc_sim_ouput <- simulateResiduals(m.aa.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.aa.g)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -3.515e+00  1.139e+00  -3.086  0.00203 **
# treatmentB   2.958e-09  7.060e-01   0.000  1.00000 

a <- Anova(m.aa.g, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 9.524  1   0.002028 **
# treatment   0.000  1   1.000000

# save anova values for table:
grid.aa.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('grid','grid'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned absent','aligned absent'))


em <- emmeans(m.aa.g, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.0289 0.032 Inf   0.00318     0.217
# B         0.0289 0.032 Inf   0.00318     0.217

grid.aa.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B             1 0.706 Inf    1   0.000  1.0000

##########################
# Create objects for viz
##########################

vis.me.g <- ggpredict(m.me.g, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.me.g$contingency <- c("ME")
vis.me.g$scale <- c("grid")

vis.ap.g <- ggpredict(m.ap.g, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.ap.g$contingency <- c("SS_y")
vis.ap.g$scale <- c("grid")

vis.dl.g <- ggpredict(m.dl.g, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.dl.g$contingency <- c("DL")
vis.dl.g$scale <- c("grid")

vis.aa.g <- ggpredict(m.aa.g, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.aa.g$contingency <- c("SS_n")
vis.aa.g$scale <- c("grid")


############################
# Create objects for table

# tab1. ANOVA

grid.anova <- rbind(grid.me.anova,grid.dl.anova,grid.ap.anova,grid.aa.anova)

# tab2. Pairs significance

grid.me.contrast <- as.data.frame(grid.me.contrast)
grid.dl.contrast <- as.data.frame(grid.dl.contrast)
grid.ap.contrast <- as.data.frame(grid.ap.contrast)
grid.aa.contrast <- as.data.frame(grid.aa.contrast)

grid.contrast <- rbind(grid.me.contrast,grid.dl.contrast,grid.ap.contrast,grid.aa.contrast)

grid.contrast$contingency <- NA
grid.contrast$contingency <- c("sink","dispersal limitation", "aligned present", "aligned absent")

grid.contrast$scale <- NA
grid.contrast$scale <- 'grid'
