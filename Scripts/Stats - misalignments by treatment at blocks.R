
#-----------------------------------------------------------------------------.
# DESCRIPTION: Fitting a glms to block level data (n = 450) for both neighbor 
# treatments, to find the estimated proportions of each (mis)alignment. Then, 
# calculate confidence intervals around these proportions, and compare proportions 
# to assess significance using emmeans
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

source("Scripts/Source - MAIN fitnessdata.R")

#--------------------------------------------------------------------------.
# POPULATION SCALE (PLOT)
#--------------------------------------------------------------------------.

########################
# Data cleaning
########################

fig2B_dat <- plotlev %>%
  select(names, grid, site, contingency, treatment, species) %>%
  distinct() %>%
  dplyr::group_by(names, treatment) %>%
  dplyr::mutate(replicates = n()) %>%
  group_by(names, treatment, contingency) %>%
  dplyr::mutate(prop = n()/replicates) %>%
  ungroup() %>%
  select(-species) %>%
  distinct() %>%
  pivot_wider(.,names_from = contingency, values_from = prop)
# all NAs = 0 because I standardized by replicate
fig2B_dat[is.na(fig2B_dat)] <- 0
fig2B_dat <- pivot_longer(fig2B_dat, cols = 6:9, names_to = "contingency", values_to = "prop")

fig2B_dat$contingency<-as.factor(fig2B_dat$contingency)
fig2B_dat$site<-as.factor(fig2B_dat$site)
fig2B_dat$grid<-as.factor(fig2B_dat$grid)
fig2B_dat$names<-as.factor(fig2B_dat$names)
fig2B_dat$treatment<-as.factor(fig2B_dat$treatment)

########################
## Fit models 
########################

################
# sinks
m.me.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'ME'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.me.p)
testDispersion(m.me.p) 
plot(fitted(m.me.p), residuals(m.me.p))
hist(residuals(m.me.p)) 
loc_sim_ouput <- simulateResiduals(m.me.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.me.p)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -0.61712    0.28623  -2.156   0.0311 *  
# treatmentB   0.32543    0.08057   4.039 5.36e-05 ***

a <- Anova(m.me.p, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept)  4.6486  1    0.03108 *  
# treatment   16.3156  1  5.362e-05 ***

# save anova values for table:
plot.me.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
           scale = c('plot','plot'),
           Df =c(a$Df[1],a$Df[2]),
           P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
           predictor = c("Intercept","Data type"),
           contingency = c('sinks','sinks'))


em <- emmeans(m.me.p, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.350 0.0652 Inf     0.235     0.486
# B         0.428 0.0698 Inf     0.299     0.566

plot.me.contrast <- pairs(em)
pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B         0.722 0.0582 Inf    1  -4.039  0.0001


################
# d. limitation
m.dl.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'DL'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.dl.p)
testDispersion(m.dl.p) 
plot(fitted(m.dl.p), residuals(m.dl.p))
hist(residuals(m.dl.p)) 
loc_sim_ouput <- simulateResiduals(m.dl.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.dl.p)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -2.0918     0.2434  -8.592  < 2e-16 ***
# treatmentB   -0.5372     0.1219  -4.407 1.05e-05 ***

a <- Anova(m.dl.p, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 73.827  1  < 2.2e-16 ***
# treatment   19.423  1  1.048e-05 ***  

# save anova values for table:
plot.dl.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                         scale = c('plot','plot'),
                         Df =c(a$Df[1],a$Df[2]),
                         P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                         predictor = c("Intercept","Data type"),
                         contingency = c('dispersal limitation','dispersal limitation'))

em <- emmeans(m.dl.p, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
#  A         0.1099 0.0238 Inf    0.0712     0.166
#  B         0.0673 0.0156 Inf    0.0424     0.105

plot.dl.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B       1.71 0.209 Inf    1   4.407  <.0001

################
# a. present
m.ap.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1| grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'SS_y'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.ap.p)
testDispersion(m.ap.p) 
plot(fitted(m.ap.p), residuals(m.ap.p))
hist(residuals(m.ap.p)) 
loc_sim_ouput <- simulateResiduals(m.ap.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # great

summary(m.ap.p)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept)  -1.5191     0.1726  -8.800  < 2e-16 ***
# treatmentB   -0.7046     0.1064  -6.622 3.53e-11 ***

a <- Anova(m.ap.p, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 77.445  1  < 2.2e-16 ***
# treatment   43.856  1  3.535e-11 ***

# save anova values for table:
plot.ap.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('plot','plot'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned present','aligned present'))

em <- emmeans(m.ap.p, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment   prob     SE  df asymp.LCL asymp.UCL
# A         0.1796 0.0254 Inf    0.1350     0.235
# B         0.0976 0.0158 Inf    0.0707     0.133

plot.ap.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B          2.02 0.215 Inf    1   6.622  <.0001

################
# a. absent
m.aa.p <- glmmTMB(prop ~ treatment + (1 | names) +  (1 | site) + (1 | grid:site), 
                  data = filter(fig2B_dat, contingency %in% 'SS_n'), 
                  weights = replicates, 
                  family = binomial()) 

# Check model fit 
testZeroInflation(m.aa.p)
testDispersion(m.aa.p) 
plot(fitted(m.aa.p), residuals(m.aa.p))
hist(residuals(m.aa.p)) 
loc_sim_ouput <- simulateResiduals(m.aa.p)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
)

summary(m.aa.p)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)  
# (Intercept) -1.01862    0.22932  -4.442 8.92e-06 ***
# treatmentB   0.35685    0.08491   4.203 2.64e-05 ***

a <- Anova(m.aa.p, type = 3)
#                 Chisq Df Pr(>Chisq)    
# (Intercept) 19.730  1  8.917e-06 ***
# treatment   17.664  1  2.636e-05 *** 

# save anova values for table:
plot.aa.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2]),3),
                            scale = c('plot','plot'),
                            Df =c(a$Df[1],a$Df[2]),
                            P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2]),3),
                            predictor = c("Intercept","Data type"),
                            contingency = c('aligned absent','aligned absent'))


em <- emmeans(m.aa.p, ~treatment, type = "response") # specify green_index_scaled values
em
# treatment  prob     SE  df asymp.LCL asymp.UCL
# A         0.265 0.0447 Inf     0.187     0.361
# B         0.340 0.0510 Inf     0.248     0.446

plot.aa.contrast <- pairs(em)
# contrast odds.ratio    SE  df null z.ratio p.value
# A / B           0.7 0.0594 Inf    1  -4.203  <.0001


##########################
# Create objects for viz
##########################

vis.me.p <- ggpredict(m.me.p, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.me.p$contingency <- c("ME")
vis.me.p$scale <- c("plot")

vis.ap.p <- ggpredict(m.ap.p, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.ap.p$contingency <- c("SS_y")
vis.ap.p$scale <- c("plot")

vis.dl.p <- ggpredict(m.dl.p, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.dl.p$contingency <- c("DL")
vis.dl.p$scale <- c("plot")

vis.aa.p <- ggpredict(m.aa.p, 
                      terms = c("treatment"), 
                      type = "fe", allow.new.levels=TRUE)
vis.aa.p$contingency <- c("SS_n")
vis.aa.p$scale <- c("plot")


############################
# Create objects for table

# tab1. ANOVA

plot.anova <- rbind(plot.me.anova,plot.dl.anova,plot.ap.anova,plot.aa.anova)

# tab2. Pairs significance

plot.me.contrast <- as.data.frame(plot.me.contrast)
plot.dl.contrast <- as.data.frame(plot.dl.contrast)
plot.ap.contrast <- as.data.frame(plot.ap.contrast)
plot.aa.contrast <- as.data.frame(plot.aa.contrast)

plot.contrast <- rbind(plot.me.contrast,plot.dl.contrast,plot.ap.contrast,plot.aa.contrast)

plot.contrast$contingency <- NA
plot.contrast$contingency <- c("sink","dispersal limitation", "aligned present", "aligned absent")

plot.contrast$scale <- NA
plot.contrast$scale <- 'plot'


