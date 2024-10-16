
##################################
# SUPPLEMENTAL MATERIAL
##################################

# DESCRIPTION:
# Supplementary material to see how effectively we reduced the cover in cleared plots in comparison to the intact paired plot. Also to see how total % cover changed across the productivity gradient.

# library(tidyverse)
# library(stringr)
# library(car)
# library(lme4)
# library(glmmTMB)
# library(visreg)
# library(patchwork)
# library(ggeffects)
# library(DHARMa)

source("Scripts/Source - MAIN fitnessdata.R")
source("Scripts/Supp figure - species specific misalignments.R")

###########################
# TOTAL ABUNDANCE DATA
###########################

# prepare plotlev abundances by name by species
plotlev$ab_cat <- as.factor(plotlev$ab_cat)
fig1_ab <- plotlev %>%
  filter(treatment %in% 'B') %>%
  select(names, species, occurrence, plot, grid, site) %>%
  mutate(presence = ifelse(occurrence == 'no', 0, 1))

# prepare species proportion of contingencies
fig1P_dat$species <- 'plaere'
fig1V_dat$species <- 'vulmic'
fig1M_dat$species <- 'miccal'
fig1B_dat$species <- 'brohor'
fig1B_dat <- select(fig1B_dat, -replicates)
spp_prop_dat <- rbind(fig1P_dat,fig1V_dat,fig1M_dat,fig1B_dat)
spp_prop_dat$species <- as.factor(spp_prop_dat$species)

# link together plotlev and proportion of contingencies for each species
review_dat <- full_join(fig1_ab, spp_prop_dat, by = c('species','names'), relationship = 'many-to-many')

# rank abundance step
# to organize x axis in species from most abundant to least:
review_dat <- na.omit(review_dat)
review_dat %>%
  dplyr::select(presence, species) %>%
  dplyr::group_by(species) %>%
  dplyr::summarize(sum_presence = sum(presence))
# 1 brohor           696
# 2 miccal          1032
# 3 plaere           836
# 4 vulmic          1048  
rank <- c('Festuca', 'Micropus', 'Plantago', 'Bromus') # most to least presences

# pull out each contingency and model:
mr1 <- glmmTMB(prop ~ species*contingency + (1|site/grid/plot), 
               data = review_dat,
               family = binomial()
               )

vis_mr1 <- ggpredict(mr1,
                     terms = c("species", "contingency"),
                     type = 'fixed'); plot(vis_mr1)

vis_mr1$x <- case_match(vis_mr1$x, #species
                     "brohor" ~ "Bromus",
                     "miccal" ~ "Micropus",
                     "plaere" ~ "Plantago",
                     "vulmic" ~ "Festuca")
vis_mr1$x <- factor(vis_mr1$x, levels = c(rank)) # by rank
vis_mr1$group <- case_match(vis_mr1$group, #contingency
                            "SS_y" ~ "Aligned present",
                            "ME" ~ "Sink",
                            "DL" ~ "Dispersal limitation",
                            "SS_n" ~ "Aligned absent")
vis_dd$group <- factor(vis_dd$group, levels = c("Aligned present", "Sink","Dispersal limitation", "Aligned absent"))

library(viridis)

review_fig <- ggplot(data = vis_mr1, mapping = aes(x = x, y = predicted, color = x)) +
  geom_point(size = 3) + 
  geom_line(color = "lightgrey", aes(group = group)) +
  geom_linerange(aes(y = predicted, ymin = conf.low, ymax = conf.high)) +
  facet_wrap(~group) +
  labs(x = "Species by rank abundance", y = "Proportion of each (mis)alignment", color = "") +
  scale_color_viridis(option = "H", direction = 1, discrete = TRUE,
                      labels = c("Bromus","Plantago", "Micropus","Festuca")) +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size = 16))

#---------------------------------------------
jpeg('Figures/fig_supp_review.jpeg', width = 6, height = 5, units = 'in', res = 600)
review_fig
dev.off()
#---------------------------------------------

# dispersal limitation
review_dl_dat




#-------------------------------

# per misalignment
# had to remove some random effects because of singular fit

mr2 <- lmer(seed ~ ab_cat*species + (1|site), data = temp_dat)
mr3 <- lmer(seed ~ poly(ab_cat,2, raw = TRUE)*species + (1|site), data = temp_dat)
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