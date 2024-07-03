
#-----------------------------------------------------------------------------
# RADAR PLOTS FOR THE CONTINGENCY PROPORTIONS CALCULATED WITH MULTINOMIALS
#-----------------------------------------------------------------------------

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
library(patchwork)

source("Scripts/Source - MAIN fitnessdata.R")
source("Scripts/Stats - misalignments by treatments by blocks.R")
source("Scripts/Stats - misalignments by treatments by grids.R")
source("Scripts/Stats - misalignments by treatments by sites.R")

# using estimated proportions for radar plots:
output2 # plot
output3 # grid
output4 # site (edit so that DL and SS_n are 0 in biotic)

# turn characters to factors
output2 <- as.data.frame(output2)
output3 <- as.data.frame(output3)
output4<- as.data.frame(output4)

output2$contingency <- as.factor(output2$contingency)
output3$contingency <- as.factor(output3$contingency)
output4$contingency <- as.factor(output4$contingency)

output2$treatment<- as.factor(output2$treatment)
output3$treatment<- as.factor(output3$treatment)
output4$treatment<- as.factor(output4$treatment)

# change names of contingencies
output2$contingency <- recode_factor(output2$contingency,
                                     'DL' = "Dispersal limited",
                                     'ME' = 'Sink population',
                                     'SS_n' = 'Aligned \n absent',
                                     'SS_y' = 'Aligned \n present')
output3$contingency <- recode_factor(output3$contingency,
                                     'DL' = "Dispersal limited",
                                     'ME' = 'Sink population',
                                     'SS_n' = 'Aligned \n absent',
                                     'SS_y' = 'Aligned \n present')
output4$contingency <- recode_factor(output4$contingency,
                                     'DL' = "Dispersal limited",
                                     'ME' = 'Sink population',
                                     'SS_n' = 'Aligned \n absent',
                                     'SS_y' = 'Aligned \n present')
# wrangle into what ggradar needs:
radar_plot <- output2 %>%
  select(treatment, prob, contingency) %>%
  pivot_wider(., names_from = contingency, values_from = prob)

radar_grid <- output3 %>%
  select(treatment, prob, contingency) %>%
  pivot_wider(., names_from = contingency, values_from = prob)

radar_site <- output4%>%
  select(treatment, prob, contingency) %>%
  pivot_wider(., names_from = contingency, values_from = prob)

# install.packages("devtools")
# devtools::install_github("ricardo-bion/ggradar", 
#                          dependencies = TRUE)
library(ggradar)

cols <- c("sienna","seagreen3")

onem2 <- ggradar(radar_plot,
                 axis.label.size = 4,
                 grid.min =0,
                 grid.max = 0.6,
                 legend.position = 'none',
                 background.circle.colour = "white",
                 axis.line.colour = "gray60",
                 gridline.min.colour = "gray60",
                 gridline.mid.colour = "gray60",
                 gridline.max.colour = "gray60",
                 group.colours = cols,
                 grid.label.size = 0,
                 group.point.size = 3) +
  theme(
    axis.text = element_text(hjust=0.2, vjust = 1),
    plot.margin = margin(0, 1, 0, 1, 'cm')) +
  coord_cartesian(clip = "off")

gridm2 <- ggradar(radar_grid,
                  axis.label.size = 4,
                  grid.min =0,
                 grid.max = 0.6,
                 legend.position = 'none',
                 background.circle.colour = "white",
                 axis.line.colour = "gray60",
                 gridline.min.colour = "gray60",
                 gridline.mid.colour = "gray60",
                 gridline.max.colour = "gray60",
                 group.colours = cols,
                 grid.label.size = 0,
                 group.point.size = 3)+
  theme(
    axis.text = element_text(hjust=0.2, vjust = 1),
    plot.margin = margin(0, 1, 0, 1, 'cm')) +
  coord_cartesian(clip = "off")
sitem2 <- ggradar(radar_site,
                 axis.label.size = 4,
                 grid.min =0,
                 grid.max = 0.6,
                 legend.text.size = 11,
                 legend.position = 'bottom',
                 background.circle.colour = "white",
                 axis.line.colour = "gray60",
                 gridline.min.colour = "gray60",
                 gridline.mid.colour = "gray60",
                 gridline.max.colour = "gray60",
                 group.colours = cols,
                 grid.label.size = 0,
                 group.point.size = 3)+
  theme(
    axis.text = element_text(hjust=0.2, vjust = 1),
    plot.margin = margin(0, 1, 0, 1, 'cm')) +
  coord_cartesian(clip = "off") # solution to labels getting cut off found here https://stackoverflow.com/questions/71790703/ggradar-how-to-increase-space-for-long-axis-labels

png("Figures/new_radar_fig2.png", units = 'in', height =8, width = 12, res = 300)
(figa + figc + fige) / (onem2 + gridm2 + sitem2)
dev.off()


