
library(fmsb)
colors_borderA=c("#70AD47")#70AD47
colors_borderB=c('#FF0000') #FF0000

colors_inA=c("#70AD47")#70AD47
colors_inB=c('#FF0000') #FF0000

## ALL SPECIES ##
# ad SE from bootstrapping code https://www.r-graph-gallery.com/142-basic-radar-chart.html

# plotlev ####
radat_plot_A <- radat_plot_all[-4,]
radat_plot_B <- radat_plot_all[-3,]

png("radat_plot_A.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_plot_A, 
           pcol=colors_borderA,
           pfcol=scales::alpha(colors_inA, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(""))
dev.off()

png("radat_plot_B.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_plot_B, 
           pcol=colors_borderB,
           pfcol=scales::alpha(colors_inB, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(""))
dev.off()

# gridlev ####
radat_grid_A <- radat_grid_all[-4,]
radat_grid_B <- radat_grid_all[-3,]

png("radat_grid_A.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_grid_A,
           pcol=colors_borderA,
           pfcol=scales::alpha(colors_inA, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(''))
dev.off()

png("radat_grid_B.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_grid_B,
           pcol=colors_borderB,
           pfcol=scales::alpha(colors_inB, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(''))
dev.off()

# sitelev ####

radat_site_A <- radat_site_all[-4,]
radat_site_B <- radat_site_all[-3,]

png("radat_site_A.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_site_A, # polygon outline color, transparency
           pcol=colors_borderA,
           pfcol=scales::alpha(colors_inA, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(''))
dev.off()

png("radat_site_B.png", height = 4, width = 4, units = "in", res=300)
radarchart(radat_site_B, # polygon outline color, transparency
           pcol=colors_borderB,
           pfcol=scales::alpha(colors_inB, 0.5) , # polygon fill color, transparency
           plwd=6 , # polygon line thickness
           plty=1,
           cglcol="grey", # color of background gris lines
           cglty=1, # background grid line type = not dashed
           axislabcol="grey", 
           caxislabels=seq(0,20,5), # ? grid
           cglwd=0.8,# ? grid
           vlcex=1, # label size
           vlabels = c(''))
dev.off()