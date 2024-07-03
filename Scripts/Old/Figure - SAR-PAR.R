
#------------------------------------------------------------------------------------
# DESCRIPTION: Figure for comparing accumulation of SAR vs PAR(s) w/ and w/out neighbors
#------------------------------------------------------------------------------------

################################
# Load packages and functions
################################

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
library(ggdark)

source("Scripts/Source - MAIN fitnessdata.R")
source("Scripts/Source - SAR manualfit.R")


# number persisting & occuring individuals at each scale by treatment (runs = 900)
# using mat_pa, mat_pb, mat_oa, mat_pb, mat_poa, mat_pob

# Post Rachel lab meeting 3/7/2022
resp_long <- pivot_longer(resp, names_to = "cat", cols = 2:6)

# Back transform area here 
#area$area <- scale(area$area, center=F) + 1 # original rescaling of area
area <- c(area$area)
area_raw <- mat_pa %>% dplyr::select(area) %>%
  arrange(area) %>%
  distinct()
area_raw <- c(area_raw$area)
new_area <- data.frame(area_raw, area)

resp_long <- left_join(resp_long, new_area, by = "area")
# end backtransform

resp_long_a <- dplyr::filter(resp_long, cat%in%c("pa","ob","poa"))
resp_long_b <- dplyr::filter(resp_long, cat%in%c("pb","ob","pob"))

mycols <- c("#88419d","cornflowerblue","midnightblue") # for pub midnightblue instead of lightgrey
#mycols <- c("#7570b3",'#d95f02',"#e7298a") # for pres
mylabs1 <- c("potential curve","diversity (SAR)","realized curve") # realized was aligned. I kind like aligned still.

# create broken axis to zoom in
#library(ggbreak) # use scale_y_cut()

#########

PARSAR_a <- 
  ggplot(data = resp_long_a, aes(x=area_raw, y=value, col = cat)) + 
  geom_line(linewidth = 1) + # par of all species
  #dark_theme_classic() +
  theme_classic() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 16),
        text = element_text(size = 16)) +
  scale_color_manual(name = "", labels = mylabs1, values = mycols) +
  labs(y="Species", x = "Area (m square)") +
  ggtitle("Without biotic interactions") +
  xlim(0,7000)+
  ylim(2,5)
  # scale_y_continuous(breaks = c(0,3:5),
  #                    labels = c(0,3:5),
  #                    limits = c(0,5)) + # close not quite
  # scale_y_cut(breaks = c(3), which = c(2,1), scales = c(0,3), space = 0.3) 
PARSAR_a

PARSAR_b <- 
  ggplot(data = na.omit(resp_long_b), aes(x=area_raw, y=value, col = cat)) + 
  geom_line(linewidth = 1) + # par of all species
  #dark_theme_classic() +
  theme_classic() +
  theme(legend.position = c(0.6, 0.15),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.text = element_text(size = 16),
        text = element_text(size = 16))+
  scale_color_manual(name = "", labels = mylabs1, values = mycols) +
  labs(y="", x = "Area (m square)") +
  ggtitle("With biotic interactions")+
  xlim(0,7000) +
  ylim(2,5)
  # scale_y_continuous(breaks = c(0,3:5),
  #                    labels = c(0,3:5),
  #                    limits = c(0,5)) + # close not quite
  # scale_y_cut(breaks = c(3), which = c(2,1), scales = c(0,3), space = 0.3) 
PARSAR_b

# -----------------------------------.
# CONFIDENCE INTERVAL FIGURE ####

levels <- factor(c("potential abiotic",
        "potential biotic",
        "diveresity (SAR)",
        "realized abiotic",
        "realized biotic"))
 ci_c$cat <- c(levels)
 ci_z$cat <- c(levels)

mycols1 <- c("#88419d", "#88419d", "cornflowerblue","midnightblue","midnightblue")

mylabs2 <- c('diversity \n (SAR)','potential \n abiotic',"potential \n biotic" ,'realized \n abiotic','realized \n biotic')

b <- ggplot(ci_c, aes(x = cat, y = mean)) +
  geom_point(col = mycols1, size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1, linewidth = 1) +
  ylim(2,4) +
  theme_classic() +
  labs(x = "", y = "Parameter value") +
  ggtitle("Intercept (c)")+
 # ylim(0,4) +
  theme(axis.text.x = element_text(angle = 310, size = 16),
        text = element_text(size = 16))  +
  scale_x_discrete(labels = mylabs2)
b

c <- ggplot(ci_z, aes(x = cat, y = mean)) +
  geom_point(col = mycols1, size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1, linewidth = 1) +
  ylim(0,1.8) +
  theme_classic() +
  labs(x = "", y = "") +
  ggtitle("Slope (z)")+
  #ylim(0,4) +
  theme(axis.text.x = element_text(angle = 310, size = 16),
        text = element_text(size = 16)) +
  scale_x_discrete(labels = mylabs2)
c

png("Figures/SAR_PAR.png", height = 12, width = 10, units = "in", res=800)
(PARSAR_a + PARSAR_b) / (b + c) 
dev.off()

# # # Biotic only parameter CIs
# 
# # mycols2 <- c( "skyblue2","slateblue","darkblue" )
#  mycols2 <- c('#d95f02',"#7570b3","#e7298a") 
#  mylabs3 <- c("SAR","PAR presence \n only","PAR presence \n & absence")   
#  
# # Biotic
# b <- ggplot(filter(ci_c,cat%in% c("PAR B","SAR","PAR occur B")), aes(x = cat, y = mean)) +
#   geom_point(col = mycols2) +
#   geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols2) +
#   ylim(2,4) +
#   theme_classic() +
#   labs(x = "", y ="Parameter value") +
#   ggtitle("Intercept (c)")+
#   theme(axis.text.x = element_blank()) +
#   scale_x_discrete(labels = mylabs3)
# 
# c <- ggplot(filter(ci_z,cat%in% c("PAR B","SAR","PAR occur B")), aes(x = cat, y = mean)) +
#   geom_point(col = mycols2) +
#   geom_errorbar(aes(ymin = lower, ymax = upper),col = mycols2) +
#   ylim(0,1.8) +
#   theme_classic() +
#   labs(x = "", y = "Parameter value") +
#   ggtitle("Slope (z)")+
#   theme(axis.text.x = element_text(angle = 320, size = 10)) +
#   scale_x_discrete(labels = mylabs3)
# 
# 
# png("Figures/sarpar_biotic_ci.png", height = 5, width = 3, units = "in", res=300)
# b / c
# dev.off()

# # PAR occur and SAR only showing sinks
# {
# mycols1 <- c("skyblue2","darkblue","darkblue")
#  mylabs2 <- c("PAR occur A","PAR occur B", 'SAR')
#  
#  b <- ggplot(filter(ci_c, cat %in% c("SAR","PAR occur A", "PAR occur B")), aes(x = cat, y = mean)) +
#    geom_point(col = mycols1) +
#    geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1) +
#    ylim(2,4) +
#    theme_classic() +
#    labs(x = "", y = "Parameter value") +
#    ggtitle("Intercept (c)")+
#    theme(axis.text.x = element_blank()) #+
#  #scale_x_discrete(labels = mylabs2)
#  
#  c <- ggplot(filter(ci_z, cat %in% c("SAR","PAR occur A", "PAR occur B")), aes(x = cat, y = mean)) +
#    geom_point(col = mycols1) +
#    geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1) +
#    ylim(0,1.8) +
#    theme_classic() +
#    labs(x = "", y = "Parameter value") +
#    ggtitle("Slope (z)")+
#    theme(axis.text.x = element_text(angle = 320, size = 10)) 
#  # 
#  png("Figures/sarpar_ci_sinks.png", height = 5, width = 3, units = "in", res=300)
#  b / c
#  dev.off()
# }
# # PAR all species and SAR only showing dispersal lim
#  {
#  mycols1 <- c("slateblue", "slateblue", "skyblue2")
#  mylabs2 <- c("PAR A","PAR B", 'SAR')
#  
#  b <- ggplot(filter(ci_c, cat %in% c("PAR A", "PAR B", "SAR")), aes(x = cat, y = mean)) +
#    geom_point(col = mycols1) +
#    geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1) +
#    ylim(2,4) +
#    theme_classic() +
#    labs(x = "", y = "Parameter value") +
#    ggtitle("Intercept (c)")+
#    theme(axis.text.x = element_blank()) #+
#  #scale_x_discrete(labels = mylabs2)
#  
#  c <- ggplot(filter(ci_z, cat %in% c("PAR A", "PAR B", "SAR")), aes(x = cat, y = mean)) +
#    geom_point(col = mycols1) +
#    geom_errorbar(aes(ymin = lower, ymax = upper), col = mycols1) +
#    ylim(0,1.8) +
#    theme_classic() +
#    labs(x = "", y = "Parameter value") +
#    ggtitle("Slope (z)")+
#    theme(axis.text.x = element_text(angle = 320, size = 10)) 
#  # 
#  png("Figures/sarpar_ci_dl.png", height = 5, width = 3, units = "in", res=300)
#  b / c
#  dev.off()
#  }
 
#--------------------------------------------------------------------------.


# What is this error? 
# Error: memory exhausted (limit reached?)
# just plot a random 50 lines in the backgroup to give an idea of variation
# or add the random effect to stat-Smooth with 95% CIs
# or use predict() to make main lines with 95% CIs IF it can take a random effect
