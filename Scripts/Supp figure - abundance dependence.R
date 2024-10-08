
#-----------------------------------------------------------------------------.
# DESCRIPTION: Supplementary figure for persistence across abundance categories
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
# library(ggeffects)
# library(emmeans)

source("Scripts/Source - MAIN fitnessdata.R")

# Abundance categories show how abundance the focal species is at a given seed production value.
# Thus it is a good measure of intraspecific density dependence

###########################
# ABUNDANCE AS CATEGORY
###########################

plotlev$ab_cat <- as.factor(plotlev$ab_cat)
temp_dat <- filter(plotlev, treatment %in% 'B')

# per species
# had to remove some random effects because of singular fit

mdd2 <- lmer(seed ~ ab_cat*species + (1|site), data = temp_dat)
mdd3 <- lmer(seed ~ poly(ab_cat,2, raw = TRUE)*species + (1|site), data = temp_dat)
AIC(mdd2,mdd3)
# df      AIC
# mdd2 18 8726.206
# mdd3 14 8742.462

# linear fit is better.

summary(mdd2)
a <- Anova(mdd2, type = 3)
# Analysis of Deviance Table (Type III Wald chisquare tests)
# 
# Response: seed
# Chisq Df Pr(>Chisq)    
# (Intercept)     2.9038  1    0.08837 .  
# ab_cat          2.7186  3    0.43707    # NOT SIGNIFICANT!
# species        21.6843  3  7.588e-05 ***
# ab_cat:species  8.2028  9    0.51385  


############# 
## STATS ##
#############

# compare category 0 to 4
temp_dat1 <- filter(temp_dat, ab_cat%in%'0' | ab_cat%in%'3')


vis_dd <- ggpredict(mdd2, 
                    terms = c("ab_cat[all]", "species"), 
                    type = "fe"); plot(vis_dd)

vis_dd$x <- case_match(vis_dd$x, 
                       "0" ~ "0",
                       "1" ~ "1-10",
                       "2" ~ "11-100",
                       "3" ~ ">101")
vis_dd$x <- factor(vis_dd$x, levels = c("0", "1-10", "11-100", ">101"))

################
# FIGURE
################
library(viridis)

length(vis_dd$x)
colortreatpred <- c(rep(c("tan4","maroon4", "turquoise1","forestgreen"), times = 4))
dodge <- position_dodge(width = 0.5)
dd2 <- ggplot() +
  # geom_jitter(temp_dat, mapping = aes(x = ab_cat, y = seed, color = species),alpha = 0.1, height = 0.1) +
  # scale_color_manual(values =c("tan4","maroon4","turquoise1","forestgreen"),
  #                    labels = c("Bromus","Plantago", "Micropus","Festuca")) + # up here so that it only deals with this layer
  geom_point(data = vis_dd, aes(x = x, y = predicted, group = group, color = group), 
             size = 4,
             position = dodge
  ) +
  geom_line(data = vis_dd, aes(x = x, y = predicted, group = group, color = group), position = dodge) +
  geom_linerange(data = vis_dd, aes(
    x = x,
    y = predicted,
    group=group,
    color = group,
    ymin = conf.low,
    ymax = conf.high),
    show.legend = F,
    linewidth = 1,
    position = dodge
  ) +
  scale_color_viridis(option = "H", direction = 1, discrete = TRUE,
                       labels = c("Bromus","Plantago", "Micropus","Festuca")) +
  # scale_color_manual(values =c("tan4","maroon4","turquoise1","forestgreen"),
  #                    labels = c("Bromus","Plantago", "Micropus","Festuca")) +
  theme_bw() +
  labs(x = "Number of individuals", y = "Seed produced per capita", color = "Species") +
 # ylim(-2,5)+
  theme(text = element_text(size = 16))
dd2

###################################################################################
jpeg('Figures/fig_supp_dd_factor.jpeg', width = 6, height = 5, units = 'in', res = 600)
dd2
dev.off()
####################################################################################

############
# Table
############

vis_dd <- as.data.frame(vis_dd)
# change species codes to names
vis_dd$group <- case_match(vis_dd$group , 
                                 "plaere" ~ "Plantago",
                                 "vulmic" ~ "Festuca",
                                 "brohor" ~ "Bromus",
                                 "miccal" ~ "Micropus")
vis_dd$group <- factor(vis_dd$group , levels = c("Bromus", "Plantago", "Micropus", "Festuca"))

tab_density_estimates <- vis_dd |>
  dplyr::mutate(conf.int =paste(round(conf.low,2), round(conf.high,2), sep = ", ")) |>
  dplyr::mutate(predicted = round(predicted,2)) |>
  dplyr::select(group, x, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
              subtitle = "Table S16. Model estimates of predicted seed production for each abundance category and for each species.")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    group = 'Species',
    x = 'Abundance category',
    predicted = 'Estimated proportion',
    conf.int = "95% CI") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))
#as_latex() # exports code but don't know how to work it right in latex

tab_density_estimates |>
  gtsave(paste0(here::here(),"/Tables/16tab_density_estimates.pdf")) 

################
# Anovas
a

d.anova <- data.frame(Chi.squared = round(c(a$Chisq[1],a$Chisq[2], a$Chisq[3], a$Chisq[4]),3),
                      Df =c(a$Df[1],a$Df[2],a$Df[3],a$Df[4]),
                      P_value =  round(c(a$`Pr(>Chisq)`[1],a$`Pr(>Chisq)`[2],a$`Pr(>Chisq)`[3],a$`Pr(>Chisq)`[4]),3),
                      predictor = c("Intercept","Abundance category","Species","Abundance category:Species")
)


all_anovas <- d.anova |>
  dplyr::select(predictor, Chi.squared, Df, P_value) |>
  gt() |>
  tab_header( title = "",
              subtitle = "Table S17. ANOVA outputs of predictors main effects and interactions. Predictors included 'abundance categories' with four levels: 0 individuals, between 1 and 10 individuals, 11-100 individuals, and over 101 individuals, and predictor species with four levels: Plantago, Bromus, Micropus, Festuca.")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    predictor = 'Predictor',
    P_value = "P-value"
  ) |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

all_anovas |>
  gtsave(paste0(here::here(),"/Tables/17tab_density_anova.pdf")) 

######################
# pairwise contrasts

em <- emmeans(mdd2, ~ab_cat | species, type = "response") # specify green_index_scaled values
d.contrast <- as.data.frame(pairs(em))
#view(d.contrast)

# change species codes to names
d.contrast$species <- case_match(d.contrast$species, 
                       "plaere" ~ "Plantago",
                       "vulmic" ~ "Festuca",
                       "brohor" ~ "Bromus",
                       "miccal" ~ "Micropus")
d.contrast$species <- factor(d.contrast$species, levels = c("Bromus", "Plantago", "Micropus", "Festuca"))

# change abundance categories to numbers
d.contrast$contrast <- case_match(d.contrast$contrast, 
                       "ab_cat0 - ab_cat1" ~ "0 / 1-10",
                       "ab_cat0 - ab_cat2" ~ "0 / 11-100",
                       "ab_cat0 - ab_cat3" ~ "0 / >101",
                       "ab_cat1 - ab_cat2" ~ "1-10 / 11-100",
                       "ab_cat1 - ab_cat3" ~ "1-10 / >101",
                       "ab_cat2 - ab_cat3" ~ "11-100 / >101")
d.contrast$contrast <- factor(d.contrast$contrast, levels = c("0 / 1-10", "0 / 11-100", "0 / >101", "1-10 / 11-100", "1-10 / >101", "11-100 / >101"))


tab_density_contrasts <- d.contrast |>
  dplyr::select(species, contrast, t.ratio, p.value) |>
  # dplyr::mutate(odds.ratio = round(odds.ratio,3)) |>
  dplyr::mutate(t.ratio = round(t.ratio,3)) |>
  dplyr::mutate(p.value = round(p.value,3)) |>
  gt() |>
  tab_header( title = "",
              subtitle = "Table S18. Contrasts between abundance categories for each species. Non significant comparisons suggest limited effects of density dependence on seed production (i.e., populations are below carrying capacity).")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    p.value = "P-value",
    contrast = 'Pairwise comparison',
    species = "Species",
    t.ratio = "t-ratio") |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))

tab_density_contrasts |>
  gtsave(paste0(here::here(),"/Tables/18tab_density_contrasts.pdf")) 



