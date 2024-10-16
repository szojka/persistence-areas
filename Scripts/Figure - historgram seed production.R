#----------------------------------------------
# RAW SEED PRODUCTION PER SPECIES
#----------------------------------------------


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
# library(sf)

source(paste0(here::here(),"/Scripts/Source - MAIN fitnessdata.R"))

seed_dat <- plotlev

# filter out germination == 0
seed_dat$status <- as.factor(seed_dat$status) # "no fruit"
levels(seed_dat$status) 

# filter for normal > 0 and no fruit and herbivory (i.e., situations where plant germinated)
seed_dat_germ1 <- seed_dat %>%
  filter(status %in% c("normal") & seed > 0)
seed_dat_germ2 <- seed_dat %>%
  filter(status %in% c("no fruit", "herbivory"))
seed_dat_germ <- rbind(seed_dat_germ1, seed_dat_germ2)   

seed_dat_germ$treatment <- case_match(seed_dat_germ$treatment, 
                               'A' ~ 'without neighbors',
                               'B' ~ 'with neighbors')

seed_dat_germ$species <- case_match(seed_dat_germ$species, 
                                 'plaere' ~ 'Plantago',
                                 'brohor' ~ 'Bromus',
                               'vulmic' ~ 'Festuca',
                               'miccal' ~ 'Micropus')

fig_seed <- ggplot(seed_dat_germ, aes(x = seed, fill = treatment, color = treatment)) +
  geom_histogram(alpha = 0.75, bins = 60) +
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c("mediumpurple1","honeydew4"), guide = FALSE) +
  labs(fill = "", x = "Seed production of germinated transplants", y = 'Count') +
  theme_bw() +
  facet_wrap(~species) +
  theme(legend.position = 'top',
        text = element_text(size = 16)) +
  geom_vline(xintercept = 2, color = 'red')
fig_seed

jpeg('Figures/supp_fig_seedprod.jpeg', width = 8.5, height = 7, units = 'in', res = 600)
fig_seed
dev.off()

#######################
# how much of the data do we cut off after 20 seeds?
tot_n <- seed_dat_germ %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
filt_n <- seed_dat_germ %>%
  filter(seed <= 20) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::select(species, freq) %>%
  distinct()
perc_p <- filt_n[,2]/tot_n[,2]

1 - perc_p  
# 1 0.13074205 Plantago
#2 0.02325581 Micropus
#3 0.02422145 Bromus
#4 0.11869436 Festuca

dat_text <- data.frame(
  label = c("2.4% of data omitted", "11.9% of data omitted", "2.3% of data omitted", "13.1% of data omitted"),
  species   = c('Bromus', 'Festuca', 'Micropus', 'Plantago'),
  treatment = 'without neighbors'
)
# p + geom_text(
#   data    = dat_text,
#   mapping = aes(x = -Inf, y = -Inf, label = label),
#   hjust   = -0.1,
#   vjust   = -1
# )

# zoomed in:
fig_seed_zoomed <- ggplot(seed_dat_germ, aes(x = seed, fill = treatment, color = treatment)) +
  geom_histogram(alpha = 0.75, bins = 60) +
  scale_fill_manual(values = c( "mediumpurple1","honeydew4")) +
  scale_color_manual(values = c("mediumpurple1","honeydew4"), guide = FALSE) +
  labs(fill = "", x = "Seed production of germinated tranplants", y = 'Count') +
  theme_bw() +
  facet_wrap(~species) +
  theme(legend.position = 'top',
        text = element_text(size = 16),
        plot.title = element_text(size = 16)) +
  geom_vline(xintercept = 2, color = 'red') +
  xlim(0, 20) + 
  geom_text(
    data    = dat_text,
    mapping = aes(x = 13, y = 80, label = label)
    #hjust   = -0.1,
    #vjust   = -1
  )
fig_seed_zoomed

jpeg('Figures/supp_fig_seedprod_zoomed.jpeg', width = 8.5, height = 7, units = 'in', res = 600)
fig_seed
dev.off()


#---------------------------------
# supporting stats
#---------------------------------

#What's the max seed production observed per species per treatment?
#I.e. it looks like the without nighbors maybe had the high seed production. So we could also mention that under realistic conditions (with neighbors) we never observed greater than ZZZ seeds produced

seed_dat_germ %>%
  dplyr::group_by(species, treatment) %>%
  dplyr::mutate(max_seed = max(seed)) %>%
  dplyr::select(species, treatment, max_seed) %>%
  distinct()
# species  treatment         max_seed
# 1 Plantago with neighbors          52
# 2 Plantago without neighbors       85*

# 3 Micropus without neighbors       35*
# 4 Micropus with neighbors          12

# 5 Bromus   with neighbors          15
# 6 Bromus   without neighbors       70*

# 7 Festuca  without neighbors       82*
# 8 Festuca  with neighbors          36
