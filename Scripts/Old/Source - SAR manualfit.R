
#-----------------------------------------------------------------------------.
# DESCRIPTION: Fit accumulation curves myself
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

# NOTE: Keep herbivory data in mat_pb and mat_ob?? Answer: see extra script that shows no qualitative difference

# I'm having asymptote problems with package mmSAR

# SAR_ALL DATAFRAME HAS NEGATIVE VALUES BECAUSE COEF DATSETS HAVE NEGATIVE Cs: problem in both without 0s post removal of herbivory and with 0s

#SOLVED:
# TROUBLESHOOT HAVING SEMI-LOG MODELS: done, made sar_exp function
# TO GET C TO BE INTERCEPT, ADD AREA = 0 TO RESPONSE DATAFRAMES: not necessary, sar intercept is put in equation right away, so new 'intercept' still reflects differences from the model coef intercepts

# Dataframes I use: mat_pa, mat_pb, mat_ob

#------------------------------------------------------------------------------.
# CALCULATE C & Z
# use lists that should be the same length as the number of runs for x(area) & s
# e.g. V11 should have one z, one c, x<- mat_pa%>%filter(run=="V11")%>% select(area)%>% distinct(), s<-mat_pa%>%filter(run == "V11")%>%select(sum)%>%distinct()

## check that the fit r2 is better with log-log or semi-log model

## mat_pa
# s <-mat_pa %>% select(sum, area, run) %>% dplyr::filter(run=='99')
# hist(s$sum)
# # semi log
# z1 <- lm(sum~log(area), data = s)
# summary(z1) # r2 40
# # log log
# z2 <- lm(log(sum+1)~log(area), data = s)
# summary(z2) # r2 in 33
# 
# ## mat_pb
# s <-mat_pb%>%dplyr::select(sum, area, run) %>% dplyr::filter(run%in%c('100'))
# hist(s)
# # semi log
# z1 <- lm(sum~log(area), data = s)
# summary(z1) # r2 is 46
# # log log
# z2 <- lm(log(sum+1)~log(area), data = s)
# summary(z2) # r2 in 38
# 
# ## mat_ob
# s <-mat_ob%>%dplyr::select(sum, area, run) %>% dplyr::filter(run%in%c('221'))
# hist(s)
# # semi log
# z1 <- lm(sum~log(area), data = s)
# summary(z1) # r2 is 1
# # log log
# z2 <- lm(log(sum)~log(area), data = s)
# summary(z2) # r2 in 0.9

# Which model is better?
#semi-log model consistently fits better

# Q Should I fit the below power_sars with the logs?
# Yes, don't remove sum == 0 observations, +1 to sum

# Q Is it biologically reasonable to remove these zeros? 
# no, it is less common to encounter 0s with rull SARs, but they would still be included if a site had no species

#------------------------------------------------------------------------------.
# FIT OVERALL MODELS ####
# with runs as random intercept and slope
# goal is to extract intercept (c) and accumulation rates (z)

# WEIRD WARNING after removing herbivory from mat_pa/pb/ob such that random slopes won't run (1+log(area)|run)
# Warning messages:
#   1: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                     Model failed to converge with max|grad| = 0.0604838 (tol = 0.002, component 1)
#   2: In checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv,  :
#                     Model is nearly unidentifiable: very large eigenvalue - Rescale variables?

# potential solutions:
# choose sampling limit in lmer formula
# rescale area to be smaller by centering and normalizing (#x$area <- rescale(x$area, to = c(0.01,1)) worked to eliminate warning #2
# for first warning use troubleshooting that Jesse helped me with
# BOlker answers; https://stackoverflow.com/questions/53034261/warning-lme4-model-failed-to-converge-with-maxgrad helped resolve warning #1

# actual solution:
#x$area <- scale(x$area, center=F) + 1 # rescale area
#mpa <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 

# do the 'fit' line as its own model for each facet with (1+log(area)|run) (make sure it's random intercept & slope)

## Fundamental PAR ####

#x <-mat_pa%>% dplyr::select(area,run,sum)
# area hist is approx normal, sum in skewed left (high concentration at high value)

## log-log model
#mpa <- lm(log(sum+1)~log(area) + run, data = x) 

## semi-log model

x <-mat_pa%>% dplyr::select(area,run,sum)
x$area <- scale(x$area, center=F) + 1 # rescale area for convergence, add 1 for logs
mpa <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 
plot(fitted(mpa), residuals(mpa))# = fine
# NOTE - when to backtransform?? solved this in 'figures_SAR-PAR.R' script
coef(mpa)
mpa.c <- coef(mpa)[1] #exp(coef(mpa)[1]) # use for log-log
mpa.z <- coef(mpa)[2]
mpa.c.se <- coef(summary(mpa))[3]
mpa.z.se <- coef(summary(mpa))[4]
mpa.n <- length(x$sum) # or 450?
rm(mpa)

## Realized PAR ####

x <-mat_pb%>% dplyr::select(area,run,sum)
x$area <- scale(x$area, center=F) + 1 # rescale area
mpb <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 
coef(mpb)
mpb.c <- coef(mpb)[1]#exp(coef(mpb)[1])
mpb.z <- coef(mpb)[2]
mpb.c.se <- coef(summary(mpb))[3]
mpb.z.se <- coef(summary(mpb))[4]
mpb.n <- length(x$sum) # of runs
rm(mpb)

## Realized SAR ####

x <-mat_ob%>% dplyr::select(area,run,sum)
x$area <- scale(x$area, center=F) + 1 # rescale area
mob <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 
coef(mob)
mob.c <- coef(mob)[1]#exp(coef(mob)[1])
mob.z <- coef(mob)[2]
mob.c.se <- coef(summary(mob))[3]
mob.z.se <- coef(summary(mob))[4]
mob.n <- length(x$sum) # or 450?
rm(mob)

# Abioitc Persistence AND Occurrence ####

x <-mat_poa%>% dplyr::select(area,run,sum)
x$area <- scale(x$area, center=F) + 1 # rescale area
mpoa <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 
plot(fitted(mpoa), residuals(mpoa))# = fine
coef(mpoa)
mpoa.c <- coef(mpoa)[1]
mpoa.z <- coef(mpoa)[2]
mpoa.c.se <- coef(summary(mpoa))[3]
mpoa.z.se <- coef(summary(mpoa))[4]
mpoa.n <- length(x$sum)
rm(mpoa)

# Bioitc Persistence AND Occurrence ####

x <-mat_pob%>% dplyr::select(area,run,sum)
x$area <- scale(x$area, center=F) + 1 # rescale area
mpob <- lmer(sum~log(area) + (1+log(area)|run), data = x, lmerControl(check.conv.grad = .makeCC("warning", tol = 0.2, relTol = NULL)), REML = F) 
plot(fitted(mpob), residuals(mpob))# = fine

coef(mpob)
mpob.c <- coef(mpob)[1]
mpob.z <- coef(mpob)[2] # why is slope above 1? ####
mpob.c.se <- coef(summary(mpob))[3]
mpob.z.se <- coef(summary(mpob))[4]
mpob.n <- length(x$sum)
rm(mpob)


#------------------------------------------------------------------------------.
# USE COEFS FROM LINEAR MODELS 
# TO FIT SAR EXP MODELS

# calculate y response for lmer models
# v1: using SAR model

# power Model: ln(S) = c + z*ln(A) for log-log 
sar_power = function(x, c, z) {
  c * (x^z) 
}

#logistic model: S = b / (c + A^-z)

#exponential model: S = c + z*ln(A) for semi-log 
# CHECK HERE - see Lauren semi-log equation ####
# use this as better with 'island' like systems, e.g. our patchy serpentine
sar_exp = function(x, c, z) {
  c + z*log(x)
}

# Replace the raw area with the min max of my scaled & transformed area
area <- mat_pa %>% dplyr::select(area) %>%
  arrange(area) %>%
  distinct()
area$area <- scale(area$area, center=F) + 1 # rescale area
x <- c(area$area)

resp <- data.frame(area = x, pa = NA, pb = NA, ob = NA, poa = NA, pob = NA)
resp$pa<- sar_exp(x, c=mpa.c,z=mpa.z)
resp$pb<- sar_exp(x, c=mpb.c,z=mpb.z)
resp$ob<- sar_exp(x, c=mob.c,z=mob.z)
resp$poa<- sar_exp(x, c=mpoa.c,z=mpoa.z)
resp$pob<- sar_exp(x, c=mpob.c,z=mpob.z)

#------------------------------------------------------------------------------.
# FIT LINEAR MODEL TO EACH RUN
# AND EXTRACT COEFS

## mat_pa
runs <- c(unique(mat_pa$run))
coef.pa <- matrix(NA, nrow=450, ncol=3)
colnames(coef.pa) <- c("c","z","run")
k <- 1
for (r in runs){
s <-mat_pa%>%filter(run == r) %>% dplyr::select(sum)
x <- mat_pa%>%filter(run == r) %>% dplyr::select(area)
z1 <- lm(s$sum~log(x$area+0.01))
summary(z1)
coef.pa[k,1] <- coef(z1)[1]
coef.pa[k,2] <- z1$coefficients[2] 
coef.pa[k,3] <- r
k <- k+1
}
coef.pa <- as.data.frame(coef.pa)
coef.pa$c <- as.numeric(coef.pa$c)
coef.pa$z <- as.numeric(coef.pa$z)

## mat_pb
coef.pb <- matrix(NA, nrow=450, ncol=3)
colnames(coef.pb) <- c("c","z","run")
k <- 1
for (r in runs){
  s <-mat_pb%>%filter(run == r)%>%dplyr::select(sum)
  x <- mat_pb%>%filter(run==r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area+0.01))
  summary(z1)
  coef.pb[k,1] <- coef(z1)[1]
  coef.pb[k,2] <-  z1$coefficients[2]
  coef.pb[k,3] <- r
  k <- k+1
}
coef.pb <- as.data.frame(coef.pb)
coef.pb$c <- as.numeric(coef.pb$c)
coef.pb$z <- as.numeric(coef.pb$z)

## mat_ob
coef.ob <- matrix(NA, nrow=450, ncol=3)
colnames(coef.ob) <- c("c","z","run")
k <- 1
for (r in runs){
  s<-mat_ob%>%filter(run == r)%>%dplyr::select(sum)
  x<- mat_ob%>%filter(run == r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area+0.01))
  summary(z1)
  coef.ob[k,1] <- coef(z1)[1]
  coef.ob[k,2] <- z1$coefficients[2]
  coef.ob[k,3] <- r
  k <- k+1
}
coef.ob <- as.data.frame(coef.ob)
coef.ob$c <- as.numeric(coef.ob$c)
coef.ob$z <- as.numeric(coef.ob$z)

## mat_poa
coef.poa <- matrix(NA, nrow=450, ncol=3)
colnames(coef.poa) <- c("c","z","run")
k <- 1
for (r in runs){
  s<-mat_poa%>%filter(run == r)%>%dplyr::select(sum)
  x<- mat_poa%>%filter(run == r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area+0.01))
  summary(z1)
  coef.poa[k,1] <- coef(z1)[1]
  coef.poa[k,2] <- z1$coefficients[2]
  coef.poa[k,3] <- r
  k <- k+1
}
coef.poa <- as.data.frame(coef.poa)
coef.poa$c <- as.numeric(coef.poa$c)
coef.poa$z <- as.numeric(coef.poa$z)

## mat_pob
coef.pob <- matrix(NA, nrow=450, ncol=3)
colnames(coef.pob) <- c("c","z","run")
k <- 1
for (r in runs){
  s<-mat_pob%>%filter(run == r)%>%dplyr::select(sum)
  x<- mat_pob%>%filter(run == r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area+0.01))
  summary(z1)
  coef.pob[k,1] <- coef(z1)[1]
  coef.pob[k,2] <- z1$coefficients[2]
  coef.pob[k,3] <- r
  k <- k+1
}
coef.pob <- as.data.frame(coef.pob)
coef.pob$c <- as.numeric(coef.pob$c)
coef.pob$z <- as.numeric(coef.pob$z)

#------------------------------------------------------------------------------.
# CALCULATE PREDICTED RESPONSES ####
## responses based on extracted c, z and predictor A

area <- mat_pa %>% dplyr::select(area) %>%
  arrange(area) %>%
  distinct()
area$area <- scale(area$area, center=F) + 1 # rescale area so its the same as what I fit the models with
x <- c(area$area)

x_length <- length(x)
runs_length <- length(runs)
y <- data.frame(rep(runs,each=x_length), rep(x, times=runs_length))
names(y)[1] <- "run"
names(y)[2] <- "area"
y$resp <- NA
y <- pivot_wider(y,names_from = run, values_from = resp)

# fill in the generic data.frame 'y':
# run in out.pa be column for 'response.pa' and initialize column y before spreading, x is another column repeated for each run
# then call run==r for c and z values in out.pa, fill y data using x 

y_pa <- y
k <- 1
for (i in 2:451){
  y_pa[,i] <- sar_exp(x, c=coef.pa[k,1],z=coef.pa[k,2])
  k <- k + 1 # counts the coef.pa row number starting from V1 for each new run
}
y_pb <- y
k <- 1
for (i in 2:451){
  y_pb[,i] <- sar_exp(x, c=coef.pb[k,1],z=coef.pb[k,2])
  k <- k+1
}
y_ob <- y
k <- 1
for (i in 2:451){
  y_ob[,i] <- sar_exp(x, c=coef.ob[k,1],z=coef.ob[k,2])
  k <- k+1
}
y_poa <- y
k <- 1
for (i in 2:451){
  y_poa[,i] <- sar_exp(x, c=coef.poa[k,1],z=coef.poa[k,2])
  k <- k+1
}
y_pob <- y
k <- 1
for (i in 2:451){
  y_pob[,i] <- sar_exp(x, c=coef.pob[k,1],z=coef.pob[k,2])
  k <- k+1
}

# meld these 5 into one dataframe: sar_all
y_ob1 <- pivot_longer(y_ob, cols = c(2:451), names_to = "run", values_to = "resp")
y_ob1$cat <- "ob"
y_pa1 <- pivot_longer(y_pa, cols = c(2:451), names_to = "run", values_to = "resp")
y_pa1$cat <- "pa"
y_pb1 <- pivot_longer(y_pb, cols = c(2:451), names_to = "run", values_to = "resp")
y_pb1$cat <- "pb"
y_poa1 <- pivot_longer(y_poa, cols = c(2:451), names_to = "run", values_to = "resp")
y_poa1$cat <- "poa"
y_pob1 <- pivot_longer(y_pob, cols = c(2:451), names_to = "run", values_to = "resp")
y_pob1$cat <- "pob"
sar_all <- rbind(y_ob1,y_pa1,y_pb1,y_poa1,y_pob1)
class(sar_all$run)
class(sar_all$cat)
sar_all$run <- as.factor(sar_all$run)
sar_all$cat <- as.factor(sar_all$cat) 
View(sar_all)
# end = this version - good for visuals

# ------------------------------------------------------------------------------.
# CONFIDENCE INTERVALS ####

# alternative way
# se and n saved from models
ci_way1 <- function(level = 0.975, n, se) qt(level,df=n-1)*se

pa_c_ci <- ci_way1(n = mpa.n, se = mpa.c.se)
pb_c_ci <- ci_way1(n = mpb.n, se = mpb.c.se)
ob_c_ci <- ci_way1(n = mob.n, se = mob.c.se)
poa_c_ci <- ci_way1(n = mpoa.n, se = mpoa.c.se)
pob_c_ci <- ci_way1(n = mpob.n, se = mpob.c.se)

pa_z_ci <- ci_way1(n = mpa.n, se = mpa.z.se)
pb_z_ci <- ci_way1(n = mpb.n, se = mpb.z.se)
ob_z_ci <- ci_way1(n = mob.n, se = mob.z.se)
poa_z_ci <- ci_way1(n = mpoa.n, se = mpoa.z.se)
pob_z_ci <- ci_way1(n = mpob.n, se = mpob.z.se)

ci_c <- data.frame(mean = c(mpa.c, mpb.c, mob.c, mpoa.c, mpob.c), upper = c(mpa.c+pa_c_ci, mpb.c+pb_c_ci, mob.c+ob_c_ci, mpoa.c+poa_c_ci, mpob.c+pob_c_ci), 
                   lower = c(mpa.c-pa_c_ci, mpb.c-pb_c_ci, mob.c-ob_c_ci, mpoa.c-poa_c_ci, mpob.c-pob_c_ci), cat = c("pa","pb","ob","poa","pob"))
ci_z <- data.frame(mean = c(mpa.z, mpb.z, mob.z, mpoa.z, mpob.z), upper = c(mpa.z+pa_z_ci, mpb.z+pb_z_ci, mob.z+ob_z_ci, mpoa.z+poa_z_ci, mpob.z+pob_z_ci), 
                   lower = c(mpa.z-pa_z_ci, mpb.z-pb_z_ci, mob.z-ob_z_ci, mpoa.z-poa_z_ci, mpob.z-pob_z_ci), cat = c("pa","pb","ob","poa","pob"))

# same qualitative results when n run = 450

# old way:

# use coef.ob, coef.pb, coef.pa
# margin for 0.025 on each end = 
# 450*0.025 # 11.25
# 450-12 # 438
# 
# # C
# 
# CI.pa.c <- coef.pa %>% dplyr::select(c) %>% arrange(c)
# CI.pa.c <- c(CI.pa.c$c)
# CI.pa.c <- CI.pa.c[12:438]
# min(CI.pa.c) # 2.053756
# max(CI.pa.c) # 3.747687
# mpa.c # 3.689778 
# 
# CI.pb.c <- coef.pb %>% dplyr::select(c) %>%arrange(c)
# CI.pb.c <- c(CI.pb.c$c)
# CI.pb.c <- CI.pb.c[12:438]
# min(CI.pb.c) # 1.424487
# max(CI.pb.c) # 3.554051
# mpb.c # 3.482092 
# 
# CI.ob.c <- coef.ob %>%  dplyr::select(c) %>% arrange(c)
# CI.ob.c <- c(CI.ob.c$c)
# CI.ob.c <- CI.ob.c[12:438]
# min(CI.ob.c) # 1.534014
# max(CI.ob.c) # 3.907241
# mob.c # 3.717467 
# 
# # Z
# 
# CI.pa.z <- coef.pa %>% dplyr::select(z) %>%arrange(z)
# CI.pa.z <- c(CI.pa.z$z)
# CI.pa.z <- CI.pa.z[12:438]
# min(CI.pa.z) # 0.02930456
# max(CI.pa.z) # 0.2330223
# mpa.z # mean = 0.428
# 
# CI.pb.z <- coef.pb %>% dplyr::select(z) %>%arrange(z)
# CI.pb.z <- c(CI.pb.z$z)
# CI.pb.z <- CI.pb.z[12:438]
# min(CI.pb.z) # 0.0536602
# max(CI.pb.z) # 0.2982476
# mpb.z # 0.7025401 
# 
# CI.ob.z <- coef.ob %>% dplyr::select(z) %>% arrange(z)
# CI.ob.z <- c(CI.ob.z$z)
# CI.ob.z <- CI.ob.z[12:438]
# min(CI.ob.z) # 0.01030525
# max(CI.ob.z) # 0.2879089
# mob.z # 0.3937409 
# 
# # dot plot with ci's instead:
# 
# ci_c <- data.frame(mean = c(mpa.c, mpb.c, mob.c), upper = c(max(CI.pa.c), max(CI.pb.c), max(CI.ob.c)), lower = c(min(CI.pa.c), min(CI.pb.c), min(CI.ob.c)), cat = c("pa","pb","ob"))
# ci_z <- data.frame(mean = c(mpa.z, mpb.z, mob.z), upper = c(max(CI.pa.z), max(CI.pb.z), max(CI.ob.z)), lower = c(min(CI.pa.z), min(CI.pb.z), min(CI.ob.z)), cat = c("pa","pb","ob"))
# 
# ggplot(ci_c, aes(x = cat, y = mean)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) +
#   ylim(0,4) +
#   theme_classic() +
#   labs(x = "", y = "")
# 
# ggplot(ci_z, aes(x = cat, y = mean)) +
#   geom_point() +
#   geom_errorbar(aes(ymin = lower, ymax = upper)) +
#   ylim(0,1) +
#   theme_classic() +
#   labs(x = "", y = "") # intervals do not touch actual value...


# ------------------------------------------------------------------------------.
# TROUBLESHOOTING

## Q why some curves are modelled to surpass 4 species 
# concluded not a problem

# 1) could it be because some runs don't start at 28 msqu?

# prob.list <- sar_all$run[sar_all$area == 28] # all of them, have to ask sooner
# mat_pa$run[mat_pa$area == 28]
# mat_pa$run[mat_pa$area == 28]
# mat_ob$run[mat_ob$area == 28]
# # for all three only one run where area = 28 "V624"
# # NOT THE PROBLEM
# 
# # 2) could be function going wrong?
# length(y_pa) #901
# length(y_pb) #901
# length(y_ob) #901
# first column is area
# NOT THE PROBLEM

#... don't think it's a problem for now

## Q Is slope constrained by upper limit in SAR facet? 
# not a stats thing, mention in discussion z could change for SAR if more species

# test <- filter(sar_all, resp > 5)
# trun <- test %>% dplyr::select(run, cat) %>% distinct()
# trun_ob <- trun%>% filter(cat == "ob")
# trun_pa <- trun%>% filter(cat == "pa") # none
# trun_pb <- trun%>% filter(cat == "pb")
# 
# check_ob <- left_join(trun_ob, coef.ob, by = 'run') #high z
# min(check_ob$z) #  0.2752599
# max(check_ob$z) # 0.4088912
# # compared to 
# min(coef.ob$z) #0.001099685
# max(coef.ob$z) #0.4088912 highest value contained in anomalies ^
# 
# check_pa <- left_join(trun_pa, coef.pa, by = 'run')
# min(check_pa$z) 
# max(check_pa$z
# # compared to 
# min(coef.pa$z) 
# max(coef.pa$z) 
# 
# check_pb <- left_join(trun_pb, coef.pb, by = 'run') #high z
# min(check_pb$z) #  0.2515631
# max(check_pb$z) # 0.4248682
# # compared to 
# min(coef.pb$z) #0.005214445
# max(coef.pb$z) #0.4248682 highest value contained in anomalies ^
