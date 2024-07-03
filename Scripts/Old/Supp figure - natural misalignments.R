
#-----------------------------------------------------------------------------.
# DESCRIPTION: SUPPLEMENTARY FIGURE BASED ON NATURAL CONDITIONS
# 1. BIOTIC MISALIGNMENTS ACROSS SCALES
# 2. BIOTIC SPECIES-SPECIFIC MISALIGNMENTS AT PLOT SCALE 
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
source("Scripts/Stats - natural misalignments (multinomial).R")


# NOTE I DID NOT APPLY THE NEW WAY OF USING COLUMN 'REPLICATES' FOR CONTROLLING FOR 0S AND NAs IN THE FOLLOWING SECTION
# figure_1_misalignments_biotic.R code to update. 

# I ALSO HAVE NOT USED DUSTY'S WAY OF MODELING CONTINGENCIES FOR THIS PART
# see ci_stats_persistence_occurrence.R to update

{
  
  # POPULATION
  
  col_treat_long <- c(fig1B_dat$colortreat)
  c <- ggplot() +
    ylim(0,1) +
    labs(x="", y = "Proportion of pops \n at 1 m squ") +
    geom_jitter(data = fig1B_dat, aes(x = contingency, y = prop), col = col_treat_long, alpha = 0.1, size = 1, width = 0.2, height = 0.1) +
    geom_point(data = summary_1B, aes(x = contingency, y = p), col = contingency_cols1, size = 3) +
    geom_errorbar(data = summary_1B, aes(x = contingency, y = p, ymin = low, ymax = high), col = contingency_cols1) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_blank ()) +
    scale_x_discrete(labels = contingency_labs)
  
  # GRID 
  
  #contingency
  fig1D_dat <- gridlev %>%
    select(grid, contingency, treatment, species) %>%
    filter(treatment %in% 'B') %>%
    distinct() %>%
    group_by(grid, contingency) %>%
    dplyr::mutate(prop = n()/total_n) %>%
    ungroup() %>%
    select(grid, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1D_dat$unique <- NA
  fig1D_dat$unique <- seq(1, length.out = length(fig1D_dat$grid))
  
  fig1D_dat1 <- fig1D_dat
  fig1D_dat1[is.na(fig1D_dat1)] <- 0
  fig1D_dat1$sum <- rowSums(fig1D_dat1[2:5])
  fig1D_dat2 <- filter(fig1D_dat1, sum != 1)
  fig1D_dat2[fig1D_dat2 == 0] <- NA
  fig1D_dat1 <- filter(fig1D_dat1, !unique %in% fig1D_dat2$unique)
  fig1D_dat <- rbind(fig1D_dat1,fig1D_dat2)
  fig1D_dat <- select(fig1D_dat, -sum)
  fig1D_dat <- pivot_longer(fig1D_dat, cols = c(2:5) ,names_to = 'contingency', values_to = 'prop')
  fig1D_dat <- na.omit(fig1D_dat)
  
  fig1D_dat$colortreat <- NA
  fig1D_dat$colortreat[fig1D_dat$contingency == "DL"]<- "slategrey"
  fig1D_dat$colortreat[fig1D_dat$contingency == "ME"]<- "violet"
  fig1D_dat$colortreat[fig1D_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1D_dat$colortreat[fig1D_dat$contingency == "SS_y"]<- "#88419d"
  
  # find CI for each contingency
  summary_1D <- fig1D_dat %>%
    dplyr::mutate(n = n()) %>% 
    # for contingencies, n should be plot n (same across groups)
    dplyr::group_by(contingency) %>%
    dplyr::mutate(p = mean(prop)) %>%
    dplyr::mutate(se = sd(prop)/sqrt(n)) %>%
    select(contingency, p, n, se) %>%
    distinct()
  summary_1D$ci <- NA
  summary_1D$ci <- ci_prop(n = summary_1D$n, p = summary_1D$p)
  summary_1D$low <- NA
  summary_1D$high <- NA
  summary_1D$low <- summary_1D$p - summary_1D$ci
  summary_1D$high <- summary_1D$p + summary_1D$ci
  summary_1D[is.na(summary_1D)] <- 0 # changing NaNs back to 0s
  summary_1D$low[summary_1D$low < 0] <- 0 # caps lower CI at 0
  
  # order so that DL, ME, SS_n, then SS_y so coresponds to colors:
  summary_1D$contingency <- as.factor(summary_1D$contingency)
  summary_1D <- summary_1D[order(summary_1D$contingency), ]
  
  col_treat_long <- c(fig1D_dat$colortreat)
  d <- ggplot() + 
    ylim(0,1) +
    labs(x="", y ="Proportion of pops \n at 25 m squ") +
    geom_jitter(data = fig1D_dat, aes(x = contingency, y = prop, col = treatment), col = col_treat_long, alpha = 0.3, size = 1.5, width = 0.2, height = 0.1) +
    geom_jitter(data = summary_1D, aes(x = contingency, y = p),position = position_dodge(width = 0.8), col = contingency_cols1, size = 3) +
    geom_errorbar(data = summary_1D, aes(x = contingency, y = p, ymin = low, ymax = high, ), col = contingency_cols1,
                  position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_blank()) +
    scale_x_discrete(labels = contingency_labs)
  
  # SITE
  
  #contingency
  fig1E_dat <- sitelev %>%
    select(site, contingency, treatment, species) %>%
    filter(treatment %in% 'B') %>%
    distinct() %>%
    group_by(site, contingency) %>%
    dplyr::mutate(prop = n()/total_n) %>%
    ungroup() %>%
    select(site, contingency, prop) %>%
    distinct() %>%
    pivot_wider(.,names_from = contingency, values_from = prop)
  
  # solve NA problem: basically proportions don't add up to 1 because of some missing data. To delineate missing data from actual 0s, I only inputed 0 in other categories if one category == 1
  fig1E_dat$unique <- NA
  fig1E_dat$unique <- seq(1, length.out = length(fig1E_dat$site))
  
  fig1E_dat1 <- fig1E_dat
  fig1E_dat1[is.na(fig1E_dat1)] <- 0
  fig1E_dat1$sum <- rowSums(fig1E_dat1[2:3])
  fig1E_dat2 <- filter(fig1E_dat1, sum != 1)
  fig1E_dat2[fig1E_dat2 == 0] <- NA
  fig1E_dat1 <- filter(fig1E_dat1, !unique %in% fig1E_dat2$unique)
  fig1E_dat <- rbind(fig1E_dat1,fig1E_dat2)
  fig1E_dat <- select(fig1E_dat, -sum)
  fig1E_dat$SS_n <- rep(0, length.out = length(fig1E_dat$site))
  fig1E_dat$DL <- rep(0, length.out = length(fig1E_dat$site))
  fig1E_dat <- pivot_longer(fig1E_dat, cols = c(2,3,5,6) ,names_to = 'contingency', values_to = 'prop')
  fig1E_dat <- na.omit(fig1E_dat)
  
  fig1E_dat$colortreat <- NA
  fig1E_dat$colortreat[fig1E_dat$contingency == "DL"]<- "slategrey"
  fig1E_dat$colortreat[fig1E_dat$contingency == "ME"]<- "violet"
  fig1E_dat$colortreat[fig1E_dat$contingency == "SS_n"]<- "#8c96c6"
  fig1E_dat$colortreat[fig1E_dat$contingency == "SS_y"]<- "#88419d"
  
  # find CI for each contingency
  summary_1E <- fig1E_dat %>%
    dplyr::mutate(n = n()) %>% 
    # for contingencies, n should be plot n (same across groups)
    dplyr::group_by(contingency) %>%
    dplyr::mutate(p = mean(prop)) %>%
    dplyr::mutate(se = sd(prop)/sqrt(n)) %>%
    select(contingency, p, n, se) %>%
    distinct()
  summary_1E$ci <- NA
  summary_1E$ci <- ci_prop(n = summary_1E$n, p = summary_1E$p)
  summary_1E$low <- NA
  summary_1E$high <- NA
  summary_1E$low <- summary_1E$p - summary_1E$ci
  summary_1E$high <- summary_1E$p + summary_1E$ci
  summary_1E[is.na(summary_1E)] <- 0 # changing NaNs back to 0s
  summary_1E$low[summary_1E$low < 0] <- 0 # caps lower CI at 0
  summary_1E$high[summary_1E$high > 1] <- 1 # caps lower CI at 1
  
  # order so that DL, ME, SS_n, then SS_y so coresponds to colors:
  summary_1E$contingency <- as.factor(summary_1E$contingency)
  summary_1E <- summary_1E[order(summary_1E$contingency), ]
  
  col_treat_long <- c(fig1E_dat$colortreat)
  e <- ggplot() + 
    ylim(0,1) +
    labs(x="", y = "Proportion of pops \n at 100 m squ") +
    geom_jitter(data = fig1E_dat, aes(x = contingency, y = prop, col = treatment), col = col_treat_long, alpha = 0.3, size = 1.5, width = 0.2, height = 0.1) +
    geom_jitter(data = summary_1E, aes(x = contingency, y = p),position = position_dodge(width = 0.8), col = contingency_cols1, size = 3) +
    geom_errorbar(data = summary_1E, aes(x = contingency, y = p, ymin = low, ymax = high, ), col = contingency_cols1,
                  position = position_dodge(width = 0.8)) +
    geom_hline(yintercept = 0.25, linetype='dashed', col = 'red') +
    theme_bw() +
    theme(legend.position = 'none',
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    scale_x_discrete(labels = contingency_labs)
}

jpeg('fig_1_scales.jpeg', width = 5, height = 5, units = 'in', res = 300)
c / d / e  
dev.off()

