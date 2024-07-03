
#-----------------------------------------------------------------------------.
# DESCRIPTION: Supplementary figure for each treatment persistence v occurrence
#-----------------------------------------------------------------------------.

# NOTES/ISSUES
# FIXME grid scale drop in persistence but rise in realized... how is this possible?

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

#-----------------------------------
# i. SHMEAR METHOD OF PERSISTENCE
source("Scripts/Stats - SAR-PAR categories across discrete scales - i. shmear.R") # loads main source too
# run first as accumulation uses same plot data.

#################
# ABIOTIC

vis1pa <- ggpredict(m.a.p, 
                      terms = c("type"), 
                      type = "fe", allow.new.levels=TRUE)
vis1pa$group <- 'plot'

vis1ga <- ggpredict(m.a.g, # FIXME predicting huge standard error, must be a mistake somewhere
                  terms = c("type"), 
                  type = "fe", allow.new.levels=TRUE)
vis1ga$group <- 'grid'

vis1sa <- ggpredict(m.a.s, 
                  terms = c("type"), 
                  type = "fe", allow.new.levels=TRUE)
vis1sa$group <- 'site'
vis1sa$conf.low[vis1sa$x %in% 'Diversity (SAR)'] <- 1

vis1a <- rbind(vis1pa,vis1ga,vis1sa)
vis1a$group <- factor(vis1a$group, c("plot","grid","site")) # order group into increasing scales
vis1a$treatment <- "Without neighbors"


#################
# BIOTIC

vis1pb <- ggpredict(m.b.p, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis1pb$group <- 'plot'

vis1gb <- ggpredict(m.b.g, # FIXME predicts huge standard error for Realized. Must be a mistake somewhere
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis1gb$group <- 'grid'

vis1sb <- ggpredict(m.b.s, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis1sb$group <- 'site'
vis1sb$conf.low[vis1sb$x %in% 'Diversity (SAR)'] <- 1

vis1b <- rbind(vis1pb,vis1gb,vis1sb)
vis1b$group <- factor(vis1b$group, c("plot","grid","site")) # order group into increasing scales
vis1b$treatment <- "With neighbors"

# remove duplicate of "Diversity (SAR)" from one of vis1b or vis1a
vis1a <- vis1a %>%
  filter(!x %in% "Diversity (SAR)")

vis1 <- rbind(vis1a, vis1b)

vis1$colour <- case_match(vis1$x, c("Potential \n without neighbors", 
                               "Potential \n with neighbors") ~ "Potential",
                     "Diversity (SAR)" ~ "Diversity (SAR)",
                     c("Realized \n with neighbors", 
                       "Diversity (SAR)","Realized \n without neighbors") ~ "Realized")

# control order that the colors appear
vis1$colour <- factor(vis1$colour,
                 levels = c("Diversity (SAR)","Potential", "Realized"))

vis1$treatment <- factor(vis1$treatment,
                 levels = c("With neighbors","Without neighbors"))

# ABIOTIC AND BIOTIC SHMEAR TOGETHER

dodge <- position_dodge(width=0.5) 
fig_par_sh <- ggplot(vis1, aes(x = group, y = predicted, color = colour, shape = treatment, linetype = treatment, group = x)) + 
  geom_point(position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(aes( color = colour, linetype = treatment, group = x), position = dodge, linewidth = 1) +
  geom_pointrange(aes(x = group, ymin = conf.low, ymax = conf.high), 
                 position = dodge) +
  theme_bw() +
  theme(legend.position = 'none',
        text = element_text(size = 16)) +
  labs(x = "Spatial scale", y = "Proportion of species", color = "", linetype = "") + # I want y to be proportion of species
  ylim(0, 1) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue", 
                                "Potential" = "orchid",
                                "Realized" = "orange"
                                )) +
  scale_linetype_manual(values = c("Without neighbors" = "dashed",
                                   "With neighbors" = "solid"), 
                        guide = "none") +
  scale_shape_manual(values = c("Without neighbors" = 2,
                                "With neighbors" = 19)) +
  guides(color=guide_legend(override.aes=list(shape=1), title = "", title.position = "left", direction = "verticle"))+ 
  guides(shape = guide_legend(title = "", title.position = "left", direction = "verticle"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 14, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20", size = 14,angle = 0, hjust = 0, face = "plain"))) +
  ggtitle("Average scaling of suitability") 
fig_par_sh

#-----------------------------------
# i. ACCUMULATION METHOD OF PERSISTENCE

source("Scripts/Stats - SAR-PAR categories across discrete scales - ii. accumulation.R") # loads main source too

#################
# ABIOTIC

vis2ap <- ggpredict(m.a.p, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis2ap$group <- 'plot'

vis2ag <- data.frame(x = c("Diversity (SAR)", "Potential \n without neighbors" ,"Realized \n without neighbors"),
                    predicted = c(1,1,1),
                    std.error = c(NA,NA,NA),
                    conf.low = c(1,1,1),
                    conf.high = c(1,1,1),
                    group = "grid")
vis2as <- data.frame(x = c("Diversity (SAR)", "Potential \n without neighbors" ,"Realized \n without neighbors"),
                    predicted = c(1,1,1),
                    std.error = c(NA,NA,NA),
                    conf.low = c(1,1,1),
                    conf.high = c(1,1,1),
                    group = "site")
vis2a <- rbind(vis2ap,vis2ag,vis2as)
vis2a$group <- factor(vis2a$group, c("plot","grid","site")) # order group into increasing scales
vis2a$treatment <- "Without neighbors"


#################
# BIOTIC

vis2bp <- ggpredict(m.b.p, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis2bp$group <- 'plot'

vis2bg <- data.frame(x = c("Diversity (SAR)", "Potential \n with neighbors" ,"Realized \n with neighbors"),
                    predicted = c(1,1,1),
                    std.error = c(NA,NA,NA),
                    conf.low = c(1,1,1),
                    conf.high = c(1,1,1),
                    group = "grid")
vis2bs <- data.frame(x = c("Diversity (SAR)", "Potential \n with neighbors" ,"Realized \n with neighbors"),
                    predicted = c(1,1,1),
                    std.error = c(NA,NA,NA),
                    conf.low = c(1,1,1),
                    conf.high = c(1,1,1),
                    group = "site")

vis2b <- rbind(vis2bp,vis2bg,vis2bs)
vis2b$group <- factor(vis2b$group, c("plot","grid","site")) # order group into increasing scales
vis2b$treatment <- "With neighbors"

# ABIOTIC AND BIOTIC ACCUMULATION TOGETHER

# remove duplicate of "Diversity (SAR)" from one of vis1b or vis1a
vis2a <- vis2a %>%
  filter(!x %in% "Diversity (SAR)")

vis2 <- rbind(vis2a, vis2b)
vis2$x <- as.character(vis2$x)
vis2$colour <- case_match(vis2$x, c("Potential \n without neighbors", 
                                    "Potential \n with neighbors") ~ "Potential",
                          "Diversity (SAR)" ~ "Diversity (SAR)",
                          c("Realized \n with neighbors", 
                            "Realized \n without neighbors") ~ "Realized")

# control order that the colors appear
vis2$colour <- factor(vis2$colour,
                      levels = c("Diversity (SAR)","Potential", "Realized"))

vis2$treatment <- factor(vis2$treatment,
                         levels = c("With neighbors","Without neighbors"))

# ABIOTIC AND BIOTIC SHMEAR TOGETHER
# FIXME not working

dodge <- position_dodge(width=0.5) 
fig_par_accum <- ggplot(vis2, aes(x = group, y = predicted, color = colour, shape = treatment, linetype = treatment, group = x)) + 
  geom_point(position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(aes(color = colour,linetype = treatment, group = x), position = dodge, linewidth = 1) +
  geom_pointrange(aes(x = group, ymin = conf.low, ymax = conf.high), 
                  position = dodge) +
  ylim(0, 1) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue", 
                                "Potential" = "orchid",
                                "Realized" = "orange"
  )) +
  scale_linetype_manual(values = c("Without neighbors" = "dashed",
                                   "With neighbors" = "solid"),
                        guide="none") +
  scale_shape_manual(values = c("Without neighbors" = 2,
                                   "With neighbors" = 19)) +
  guides(color=guide_legend(override.aes=list(shape=1), title = "", title.position = "left", direction = "verticle"))+ 
  guides(shape = guide_legend(title = "", title.position = "left", direction = "verticle"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20",size = 16, angle = 0, hjust = 0, face = "plain"))) +
  labs(x = "Spatial scale", y = "", color = "", linetype = "", title = "Accumulated scaling of suitability") +
  ggtitle("Accumulated scaling of suitability") +
  theme_bw() +
  theme(legend.position = c(0.7,0.4),
        text = element_text(size = 16)) 
fig_par_accum 

#################
# Bring together

png("Figures/discrete_SAR-PAR.png", height = 8, width = 14, units = "in", res=800)
fig_par_sh + fig_par_accum 
dev.off()


################################################
# SIMULATE CONCEPTUAL COMPONENT

t <- 18*2 # 2 scaling types * 3 div cats * 3 scales * 2 scenarios

# In real data scales = x, predicted = predicted, group = div_type
scenario <- c(rep("scenario i", times = t/2), rep("scenario ii", times = t/2))
scale_type <- c(rep(c(rep("Accumulated", times = t/4), rep("Averaged", times = t/4)), times = 2))
scales <- c(rep(c(rep("plot", times = t/12),rep("grid", times = t/12),rep("site", times = t/12)),times = 4))
div_type <- c(rep(c("Diversity (SAR)","Realized","Potential"), times = t/3))
predicted <- c(12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
              12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
              12/16, 6/16, 7/16, 1, 1, 1, 1, 1, 1,
              12/16, 6/16, 7/16, 0, 0, 0, 0, 0, 0)

# a accumulated plot diverity = 12/16
# a accumulated plot realized = 6/16
# a accumulated plot potential = 7/16

# a accumulated grid diverity = 1/1
# a accumulated grid realized = 1/1
# a accumulated grid potential = 1/1

# a accumulated site diverity = 1/1
# a accumulated site realized = 1/1
# a accumulated site potential = 1/1

# a averaged plot diverity = 12/16
# a averaged plot realized =6/16
# a averaged plot potential = 7/16

# a averaged grid diverity = 1/1
# a averaged grid realized = 1/1
# a averaged grid potential = 1/1

# a averaged site diverity = 1/1
# a averaged site realized = 1/1
# a averaged site potential = 1/1

#-----------

# b accumulated plot diverity = 12/16
# b accumulated plot realized = 6/16
# b accumulated plot potential = 7/16

# b accumulated grid diverity = 1
# b accumulated grid realized = 1
# b accumulated grid potential = 1

# b accumulated site diverity = 1
# b accumulated site realized = 1
# b accumulated site potential = 1

# b averaged plot diverity = 12/16
# b averaged plot realized =6/16
# b averaged plot potential = 7/16

# b averaged grid diverity = 0
# b averaged grid realized =0
# b averaged grid potential =0

# b averaged site diverity = 0
# b averaged site realized =0
# b averaged site potential =0

dat <- data.frame(scenario, scale_type,scales,div_type,predicted)

dat$scales <- factor(dat$scales, c("plot","grid","site"))

ggplot(dat, aes(x = scales, y = predicted, color = div_type, group = div_type)) + 
  geom_point(position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(aes(color = div_type, group = div_type), position = dodge, linewidth = 1) +
  ylim(0, 1) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue", 
                                "Potential" = "orchid",
                                "Realized" = "orange"
  )) +
  guides(color=guide_legend(override.aes=list(shape=1), title = "", title.position = "left", direction = "horizontal"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20",size = 16, angle = 0, hjust = 0, face = "plain"))) +
  labs(x = "Spatial scale", y = "Proportion of populations", color = "", linetype = "", title = "") +
  facet_grid(cols = vars(scale_type), rows = vars(scenario)) +
  theme_bw() +
  theme(legend.position = "top",
        text = element_text(size = 16)) 


vis1$scenario <- "data"
vis2$scenario <- "data"

vis2$scale_type <- "accumulated"
vis1$scale_type <- "averaged"

vis <- rbind(vis1, vis2)

dodge <- position_dodge(width=0.5) 
fig_par <- ggplot(vis, aes(x = group, y = predicted, color = colour, shape = treatment, linetype = treatment, group = x)) + 
  geom_point(position = dodge,inherit.aes = TRUE, size = 3) +
  geom_line(aes(color = colour,linetype = treatment, group = x), position = dodge, linewidth = 1) +
  geom_pointrange(aes(x = group, ymin = conf.low, ymax = conf.high), 
                  position = dodge) +
  ylim(0, 1) +
  scale_color_manual(values = c("Diversity (SAR)" = "blue", 
                                "Potential" = "orchid",
                                "Realized" = "orange"
  )) +
  scale_linetype_manual(values = c("Without neighbors" = "dashed",
                                   "With neighbors" = "solid"),
                        guide="none") +
  scale_shape_manual(values = c("Without neighbors" = 2,
                                "With neighbors" = 19)) +
  guides(color=guide_legend(override.aes=list(shape=1), title = "", title.position = "left", direction = "verticle"))+ 
  guides(shape = guide_legend(title = "", title.position = "left", direction = "verticle"))+ 
  guides(theme(legend.title = element_text(color = "black", size = 16, angle = 0, hjust = .5, face = "plain"),
               legend.text=element_text(color = "grey20",size = 16, angle = 0, hjust = 0, face = "plain"))) +
  labs(x = "Spatial scale", y = "", color = "", linetype = "", title = "Accumulated scaling of suitability") +
  theme_bw() +
  facet_grid(cols = vars(scale_type), rows = vars(scenario)) +
  theme(legend.position = c(0.7,0.4),
        text = element_text(size = 16)) 
fig_par




