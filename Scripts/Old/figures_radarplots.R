# Fig 1.2 ####
# each contingency is a node, treatment are plots (pool species), contingency y or n are min and max spread contingencies into columns  


#---------------------------------------------------------------
# Using ggradar instead:
# install.packages("devtools")
# devtools::install_github("ricardo-bion/ggradar")
library(ggradar)

# prep data into single column with contingencies (groups), and multiple columns with values

# ALL SPECIES ####

cols <- c("sienna","seagreen")

## plot level

radat_plot_all <- plotlev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n)

onem2 <- ggradar(radat_plot_all,
        grid.min = 148,
        grid.max = 694,
        legend.position = 'none',
        background.circle.colour = "white",
        axis.line.colour = "gray60",
        gridline.min.colour = "gray60",
        gridline.mid.colour = "gray60",
        gridline.max.colour = "gray60",
        group.colours = cols,
        grid.label.size = 0)

radat_grid_all <- gridlev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct()

twenty5m2 <- ggradar(radat_grid_all,
        grid.min = 1,
        grid.max = 56,
        legend.position = 'none',
        background.circle.colour = "white",
        axis.line.colour = "gray60",
        gridline.min.colour = "gray60",
        gridline.mid.colour = "gray60",
        gridline.max.colour = "gray60",
        group.colours = cols,
        grid.label.size = 0)


radat_site_all <- sitelev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  dplyr::distinct()
radat_site_all$DL <- c(0,0)
radat_site_all$SS_n <- c(0,0)
radat_site_all = radat_site_all[,c(1,4,2,5,3)] # rearrrange column order so radar the same

one00m2 <- ggradar(radat_site_all,
        grid.min = 0,
        grid.max = 21,
        legend.position = "bottom",
        background.circle.colour = "white",
        axis.line.colour = "gray60",
        gridline.min.colour = "gray60",
        gridline.mid.colour = "gray60",
        gridline.max.colour = "gray60",
        group.colours = cols,
        grid.label.size = 0)  # Affects the grid annotations (0%, 50%, etc.)

library(patchwork)

onem2 +  twenty5m2 + one00m2

# BY SPECIES ####



#---------------------------------------------------------------
#source("Scripts/radar_data.R")

library(fmsb)
colors_border=c( rgb(0.8,0.2,0.5,0.9), rgb(0.2,0.5,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.8,0.2,0.5,0.9), rgb(0.2,0.5,0.5,0.9) , rgb(0.7,0.5,0.1,0.4) )

## ALL SPECIES ##
# ad SE from bootstrapping code https://www.r-graph-gallery.com/142-basic-radar-chart.html

#plotlev
radarchart(radat_plot_all, 
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))
# asthetic help https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html

# gridlev
radarchart(radat_grid_all,
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))


# sitelev
radarchart(radat_site_all, # polygon outline color, transparency
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))


# legend
# legend(x=1.2, y=1.2, 
#        legend = rownames(radat_site_all[-c(1,2),]), 
#        bty = "n", 
#        pch=20 , 
#        col=colors_in , 
#        text.col = "grey", 
#        cex=1.2, 
#        pt.cex=3)

## BY SPECIES ##

#### miccal ####

#plotlev
radarchart(radat_plot_miccal, 
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))
# asthetic help https://www.r-graph-gallery.com/143-spider-chart-with-saveral-individuals.html

# gridlev
radarchart(radat_grid_miccal,
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))


# sitelev
radarchart(radat_site_miccal, # polygon outline color, transparency
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

#### plaere ####

#plotlev
radarchart(radat_plot_plaere, 
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))
# gridlev
radarchart(radat_grid_plaere,
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

# sitelev
radarchart(radat_site_plaere, # polygon outline color, transparency
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

#### vulmic ####

#plotlev
radarchart(radat_plot_vulmic, 
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

# gridlev
radarchart(radat_grid_vulmic,
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

# sitelev
radarchart(radat_site_vulmic, # polygon outline color, transparency
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

#### brohor ####

#plotlev
radarchart(radat_plot_brohor, 
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

# gridlev
radarchart(radat_grid_brohor,
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))

# sitelev
radarchart(radat_site_brohor, # polygon outline color, transparency
           pcol=colors_border,
           pfcol=scales::alpha(colors_in, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c('unmatched (DL)', 'unmatched \n (sink)', 'matched (-)','matched \n (+)'))
