## Paper 1 - spatial scaling

## contains figs 
# 1.1 (map of persistence vs occupancy),
source("Scales project/Scripts/Source_fitness data.R") # chain source: this script sources main

# Main pie charts for map
# summary stats - n of each type
pie_biotic <- filter(sitelev, treatment%in%'B')
pie_biotic$group <- NA
pie_biotic$contingency <- as.factor(pie_biotic$contingency)
pie_biotic$group[pie_biotic$contingency == c("ME")] <- 'unmatched'
pie_biotic$group[pie_biotic$contingency == c("SS_y")] <- 'matched'

pie_biotic <- table(pie_biotic$group, pie_biotic$site)
pie_biotic <- as.data.frame(pie_biotic)

png("pie_biotic_fig.png", height = 5, width = 8, "in", res = 300)
ggplot(pie_biotic, aes(x = "", y = Freq, fill = Var1)) +
  geom_col(color = "black") +
  #geom_text(aes(label = Freq),position = position_fill(vjust = 1)) +
  coord_polar(theta = "y") +
  scale_fill_brewer() +
  theme_void() + facet_wrap(~Var2)#+
dev.off() 



# Summary stats ####
# see df.scales for n
mary <- plotlev %>% filter(treatment %in% "A")
table(mary$contingency, mary$species)
mary <- table(mary$contingency)
mary <- as.data.frame(mary)
mary$group <- NA
mary$group[mary$Var1 == c("DL","ME")] <- 'unmatched'
mary$group[mary$Var1 == c("SS_n","SS_y")] <- 'matched'

# brohor (tot n = 337) % plots that are sinks = 23.1454. sources = 20.178%
# miccal (tot = 355) sinks = 60% sources = 3.0986%
# plaere (tot = 334) sinks = 24.85%, sources = 27.545%
# vulmic (tot = 370) sinks = 31.62%, sources = 32.973%

beth <- plotlev %>% filter(treatment %in% "B")
table(beth$contingency, beth$species)
beth <- table(beth$contingency)
beth <- as.data.frame(beth)
beth$group <- NA
beth$group[beth$Var1 == c("DL","ME")] <- 'unmatched'
beth$group[beth$Var1 == c("SS_n","SS_y")] <- 'matched'
# brohor (tot = 281) plots that are sinks = 38.79%. sources = 16.014%
# miccal (tot = 408) sinks = 50%. sources = 1.96 %
# plaere (tot = 397) sinks = 26.70 %. sources = 14.609%
# vulmic (tot = 404) sinks = 34.158%. sources = 19.059%

#-----------------------------------------------------------------------end sum stats

# WHY is line 31 here?? ####
names(utm)[1] <- "names"

#-----------------------------------------------------------------------.
#### Fig 1.1 ####

# map (simplified area somehow… visreg?) with blue->red heat gradient. 
# facet 1 = seed production, facet 2 = occurrence 

# To visualise this use a scatter plot where you have X coordinates in one column, y coordinates in another column, yes no 2 persistence, and yes no to occurrence,this creates fine scale mapping.Summary statistics as in means and variance around each physical plot battleship point could be combined to create onegred visualizacion heat map

# join with nkey (from Distance_matrix_UTMs)
nkey <- dplyr::select(key, -ab_cat,-plot,-image)
nkey$names1<-paste(nkey$site, nkey$grid, sep="")
num <- rep(c(1:25), each = 8)
nkey$names <- paste(num, nkey$names1, sep="_")
nkey <- dplyr::select(nkey, -names1)

# plotlev has yes / no data
ttemp <- dplyr::select(plotlev, tag, species, treatment, grid, site, contingency)
mad1 <- left_join(ttemp, nkey, by = c("tag", "species", "treatment", "grid", "site"))
mad <- left_join(mad1, utm, by = "names")
mad <- dplyr::select(mad, -tag,-grid,-site)
mad <- separate(mad, col = "names", into = c("map_x","name"), sep = "_")

# create x and y axis to overlay all grids
mad$map_y <- NA

mad$map_y[mad$map_x%in%c(1:5)] <- 5
mad$map_y[mad$map_x%in%c(6:10)] <- 4
mad$map_y[mad$map_x%in%c(11:15)] <- 3
mad$map_y[mad$map_x%in%c(16:20)] <- 2
mad$map_y[mad$map_x%in%c(21:25)] <- 1

mad$map_x[mad$map_x%in%c(1,10,11,20,21)] <- 1
mad$map_x[mad$map_x%in%c(2,9,12,19,22)] <- 2
mad$map_x[mad$map_x%in%c(3,8,13,18,23)] <- 3
mad$map_x[mad$map_x%in%c(4,7,14,17,24)] <- 4
mad$map_x[mad$map_x%in%c(5,6,15,16,25)] <- 5

rm(ttemp, mad1)

# plot
ggplot(mad, aes(x = map_x, y=map_y, color = contingency)) +
  geom_jitter(alpha = 0.5)+
  theme_classic() + 
  facet_wrap(~species)

#-----------------------------------------------------------------------.

#-----------------------------------------------------------------------.

#----------------------------------------------------------------.