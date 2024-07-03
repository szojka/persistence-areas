## SOURCE FITNESS DATA & SAR MANUAL FIT WITH HERBIVORY DATA

# NOTES:
# MAT_OB ONLY LENGTH TO HAVE CHANGED WITH INCLUSION OF HERBIVORY
# (1+log(area)|run) STILL WON'T CONVERGE (SAME TWO WARNINGS MPA)
# FIGURE IS SIMILAR TO OG ONE IN DRAFT, MEANING HERBIVORY CHANGES MODEL FIT FOR REALIZED DISTRIBUTION
# YET HERBIVORY ONLY KEEPS 4 OBSERVATIONS OF ZERO SEEDS BUT 2 OBS WITH POSITIVE ABUNDANCES (EXPLANS CHANGE TO REALIZED DISTRIBUTION)
# FOUND ISSUE IN REMOVING TOO MANY ZEROS FROM SEEDS_UNCLEAN, RERUNING EVERYTHING WITH HERBIVORY, AS THIS IS A MEANINGFUL ZERO LIKE SEEDS LANDING ON ROCKS

# SOURCE ####
#### Load Packages
library(readr)
library(tidyverse)
library(raster)
library(truncnorm)
library(stringr)
library(BiodiversityR)
library(car)
library(lme4)
library(glmmTMB)
library(visreg)
library(Rmisc)
library(patchwork)

seeds_unclean <- read_csv("Master data file/Transplant_master data_Oct2020.csv")

seeds_unclean$status <- as.factor(as.character(seeds_unclean$status))
levels(seeds_unclean$status) # four levels

#### Clean dataset
seeds_unclean$status <- as.character(as.factor(seeds_unclean$status)) # make status = NA into factor
seeds_unclean$status[is.na(seeds_unclean$status)] <- c("normal")
#seeds_unclean$status <- as.factor(seeds_unclean$status)
seeds_unclean <- filter(seeds_unclean, !status%in%c('missing'))
seeds_unclean <- na.omit(seeds_unclean)

#### Summarize data
# summarize by scale to get means, n, and cis grouped by treatments
seeds_clean <- ungroup(seeds_unclean)
sum_grid1 <- summarySE(seeds_clean, measurevar="seed", groupvars=c("grid","species","treatment"))
sum_site1 <- summarySE(seeds_clean, measurevar="seed", groupvars=c("site","species","treatment"))
sum_region <- summarySE(seeds_clean, measurevar="seed", groupvars=c("species","treatment"))

sum_grid1 <- dplyr::select(sum_grid1, -grid) # remove named columns
sum_site1 <- dplyr::select(sum_site1, -site)

# because there are many grid, many site, observation, no summarize as grid as the observations
sum_grid <- summarySE(sum_grid1, measurevar="seed", groupvars=c("species","treatment"))
sum_site <- summarySE(sum_site1, measurevar="seed", groupvars=c("species","treatment"))

sum_grid$scale <- 25 # add scale column to each
sum_site$scale <- 75
sum_region$scale <- 450 #(75*6)

sum_scales <- rbind(sum_grid, sum_site, sum_region)
sum_scales$scale <- factor(sum_scales$scale, c(25,75,450)) # reorder factor levels
sum_scales$scale <- as.character(as.factor(sum_scales$scale)) # make scale continuous
sum_scales$scale <- as.numeric(as.character(sum_scales$scale)) # make scale continuous

rm(seeds_unclean, sum_grid1,sum_grid,sum_site, sum_region, sum_site1)

#### Add scale column to seeds_clean
# for Figure 2 script
# I think they will have to be different datasets
seeds_ind <- dplyr::select(seeds_clean, -grid, -site)
seeds_grid <- dplyr::select(seeds_clean, -site, -tag)
seeds_site <- dplyr::select(seeds_clean, -grid, -tag)
seeds_region <- dplyr::select(seeds_clean, -grid,-site,-tag)

#### NDVI ab
# note it was collected using point locations. Might be able to get averages from surrounding grids as well
ndvi_summary <- read_csv("McLaughlin site info/ndvi_szojka_buffer25m2_06-03-2016.csv")
ndvi_mean <- dplyr::select(ndvi_summary, grid, mean_ndvi)
ndvi_seeds <- left_join(seeds_clean, ndvi_mean, c("grid")) # join with seeds_clean

# abundance to ndvi_seeds (temporarily ndvi_ab so that I don't loose seeds data)

# abundance categories = 0,1,2,3
# key = 0, 1-10, 11-100, 101-1000

# preliminary miccal and plaere only
photos1 <- read_csv("Data/Abundance & occurrence 2019_final.csv")
# gather plaere and miccal to species
photos <- gather(photos1, species, ab_cat, 6:9) # abundance categories
# join tags to photos using:
tag_labels <- read_csv("McLaughlin site info/Tag labels.csv")
key <- left_join(tag_labels, photos, key = c("species", "treatment", "grid", "site", "plot"))
key <- na.omit(key) # remove NAs (brohor and vulmic for now)
ndvi_ab1 <- left_join(ndvi_seeds, key, by =c("tag","species","treatment","grid","site"))
ndvi_ab1 <- filter(na.omit(ndvi_ab1))

# there is a single date in the ab_cat -> fix
ndvi_ab1$ab_cat <- as.numeric(ndvi_ab1$ab_cat)
# observation 1262 is culprit -> NA = 2
ndvi_ab1$ab_cat[is.na(ndvi_ab1$ab_cat)] <- 2

ndvi_ab1$ab_cat <- as.numeric(as.character(ndvi_ab1$ab_cat))

plots <- levels(factor(ndvi_ab1$plot))
grids <- as.numeric(levels(factor(ndvi_ab1$grid)))
spp <- levels(factor(ndvi_ab1$species))

ndvi_ab1$occ <- ndvi_ab1$ab_cat

ndvi_ab1 %>% filter(species == "plaere") %>% na.omit(ab_cat) -> ndvi_ab2

for(j in grids) {
  plots <- levels(factor(ndvi_ab2$plot[ndvi_ab2$grid == j]))
  for(i in plots) {
    if(length(ndvi_ab2$ab_cat[ndvi_ab2$treatment %in% c("A","B") & ndvi_ab2$plot == i & ndvi_ab2$grid == j]) == 2) {
      if(ndvi_ab2$ab_cat[ndvi_ab2$treatment == "A" & ndvi_ab2$plot == i & ndvi_ab2$grid == j] <= ndvi_ab2$ab_cat[ndvi_ab2$treatment == "B" & ndvi_ab2$plot == i & ndvi_ab2$grid == j]){
        
        ndvi_ab2$occ[ndvi_ab2$treatment == "A" & ndvi_ab2$plot == i & ndvi_ab2$grid == j] <- ndvi_ab2$occ[ndvi_ab2$treatment == "B" & ndvi_ab2$plot == i & ndvi_ab2$grid == j]
      }
    }
  }
}

# run loop again for each species then rbind
ndvi_ab1 %>% filter(species == "miccal") %>% na.omit(ab_cat) -> ndvi_ab3

for(j in grids) {
  plots <- levels(factor(ndvi_ab3$plot[ndvi_ab3$grid == j]))
  for(i in plots) {
    if(length(ndvi_ab3$ab_cat[ndvi_ab3$treatment %in% c("A","B") & ndvi_ab3$plot == i & ndvi_ab3$grid == j]) == 2) {
      if(ndvi_ab3$ab_cat[ndvi_ab3$treatment == "A" & ndvi_ab3$plot == i & ndvi_ab3$grid == j] <= ndvi_ab3$ab_cat[ndvi_ab3$treatment == "B" & ndvi_ab3$plot == i & ndvi_ab3$grid == j]){
        
        ndvi_ab3$occ[ndvi_ab3$treatment == "A" & ndvi_ab3$plot == i & ndvi_ab3$grid == j] <- ndvi_ab3$occ[ndvi_ab3$treatment == "B" & ndvi_ab3$plot == i & ndvi_ab3$grid == j]
      }
    }
  }
}


ndvi_ab1 %>% filter(species == "brohor") %>% na.omit(ab_cat) -> ndvi_ab4

for(j in grids) {
  plots <- levels(factor(ndvi_ab4$plot[ndvi_ab4$grid == j]))
  for(i in plots) {
    if(length(ndvi_ab4$ab_cat[ndvi_ab4$treatment %in% c("A","B") & ndvi_ab4$plot == i & ndvi_ab4$grid == j]) == 2) {
      if(ndvi_ab4$ab_cat[ndvi_ab4$treatment == "A" & ndvi_ab4$plot == i & ndvi_ab4$grid == j] <= ndvi_ab4$ab_cat[ndvi_ab4$treatment == "B" & ndvi_ab4$plot == i & ndvi_ab4$grid == j]){
        
        ndvi_ab4$occ[ndvi_ab4$treatment == "A" & ndvi_ab4$plot == i & ndvi_ab4$grid == j] <- ndvi_ab4$occ[ndvi_ab4$treatment == "B" & ndvi_ab4$plot == i & ndvi_ab4$grid == j]
      }
    }
  }
}

ndvi_ab1 %>% filter(species == "vulmic") %>% na.omit(ab_cat) -> ndvi_ab5

for(j in grids) {
  plots <- levels(factor(ndvi_ab5$plot[ndvi_ab5$grid == j]))
  for(i in plots) {
    if(length(ndvi_ab5$ab_cat[ndvi_ab5$treatment %in% c("A","B") & ndvi_ab5$plot == i & ndvi_ab5$grid == j]) == 2) {
      if(ndvi_ab5$ab_cat[ndvi_ab5$treatment == "A" & ndvi_ab5$plot == i & ndvi_ab5$grid == j] <= ndvi_ab5$ab_cat[ndvi_ab5$treatment == "B" & ndvi_ab5$plot == i & ndvi_ab5$grid == j]){
        
        ndvi_ab5$occ[ndvi_ab5$treatment == "A" & ndvi_ab5$plot == i & ndvi_ab5$grid == j] <- ndvi_ab5$occ[ndvi_ab5$treatment == "B" & ndvi_ab5$plot == i & ndvi_ab5$grid == j]
      }
    }
  }
}


ndvi_ab <- rbind(ndvi_ab2, ndvi_ab3, ndvi_ab4, ndvi_ab5)
names(ndvi_ab)[12] <- "occurrence" 

rm(photos1,ndvi_ab1,ndvi_ab2, ndvi_ab3, ndvi_ab4, ndvi_ab5, ndvi_seeds)

# scale occurrence from 0-1 (numeric)
library(scales)
ndvi_ab$occurrence <- rescale(ndvi_ab$occurrence, c(0,1))

#### UTM & dat

utm <- read_csv("Data/UTM_szojkalocs_id.csv")
utm2 <- dplyr::select(utm, -name) 
utm2 <- as.matrix(utm2)

# Random Location Sampling

max(utm$utm_x) #557021.5
min(utm$utm_x) #547533.3

max(utm$utm_y) #4302411
min(utm$utm_y) #4297447

{set.seed(123)
  nrun <- 900 # the number of runs
  utm_x <- rtruncnorm(n=nrun, a=547533.3, b=557021.5, mean=552277.4, sd=4000)
  utm_y <- rtruncnorm(n=nrun, a=4297447, b=4302411, mean=4299929, sd=2000)
}

pts <- data.frame(utm_x,utm_y) # dataset of 99 random points in my landscape
pts <- as.matrix(pts)

# find distances between each focal point in "pts" and my plot points in "utm"\
x <- pointDistance(pts, utm2, lonlat= FALSE, allpairs = TRUE)

# distance dataset
names <- c(utm$name)# create vector of names from utm

dist.dat = as.data.frame(t(x))

# Should retain order such that unique key (names) is maintained for each group

# make a key to connect 'names' to 'seeds_clean' (Step 2 & 3)
nkey <- dplyr::select(key, -ab_cat,-plot,-image)
nkey$names1<-paste(nkey$site, nkey$grid, sep="")
num <- rep(c(1:25), each = 8)
nkey$names <- paste(num, nkey$names1, sep="_")
nkey <- dplyr::select(nkey, -names1)

# change seed column to persistence
pers <- dplyr::select(seeds_clean, -status)
#pers <- pivot_wider(pers, names_from = "species", values_from = "seed")
name_key <- left_join(pers, nkey, by = c("tag","species","treatment","grid","site"))
name_key <- na.omit(name_key)

# temporary filter to make persistence column
# for each species, in each 'unique' value, sumarize into new column
name_key$unique <- paste(name_key$names, name_key$treatment, sep = "_")

p_sum <- name_key %>%
  dplyr::select(species, seed, unique)

p_sum$seed[p_sum$seed > 0] <- 1

plot_list <- p_sum %>% 
  group_by(unique) %>% 
  add_count() %>% 
  filter(n == 4)

uni <- p_sum %>% 
  filter(unique %in% plot_list$unique) %>%
  distinct() %>%
  pivot_wider(names_from = "species", values_from = "seed") %>%
  mutate(persist = rowSums(.[2:5]))

uni_key1 <- dplyr::select(uni, unique, persist)
uni_key <- left_join(uni_key1, name_key, by = "unique")
uni_key <- uni_key %>%
  dplyr::select(-unique, -seed, -tag, -species)
rm(uni, uni_key1, p_sum, plot_list, pers)

# collapse dataframe (dist.dat) into two columns where column names are third column indicating value 
dat1 <- pivot_longer(dist.dat, cols = c(1:nrun), names_to = "run")
dat1$names <- rep(names, each=nrun)
names(dat1)[2] <- c("area")
dat <- left_join(dat1, uni_key, "names")
dat <- distinct(dat)

# checking all is well with this dataset
str(dat)
dat$run <- as.factor(as.character(dat$run))
levels(dat$run) # good: 450 exist
dat$treatment <- as.factor(as.character(dat$treatment))
levels(dat$treatment) # good
dat <- na.omit(dat)

# problem with statistical modelling as area and percent of persisting species are on way difference scales
# change area to km (/1000)
dat$area <- dat$area/1000

# scale category for joy divisions
dat$scale <- cut(dat$area, breaks = 15, labels = c(0:14), include.lowest = T)
dat$scale <- as.numeric(as.factor(dat$scale)) # do this so my loop below will run
table(dat$scale) # frequency of obs w/in each level is well distributed

dat$run <- as.character(as.factor(dat$run))

#### plotlev

# abundance categories = 0,1,2,3
# ab_cat = 0, 1-10, 11-100, 101-1000

# join seeds_clean to photo abundances
plotlev1 <- left_join(key, seeds_clean, key = c("tag", "species", "treatment", "grid", "site"))
plotlev1 <- dplyr::filter(plotlev1, !status%in%NA) # NAs represent observations I do not have

# add occurrence vs persistence columns
#plotlev1 <- plotlev1 %>%
#mutate(occurrence = ab_cat) %>%
#mutate(persistence = seed)
#plotlev <- plotlev1 # allows plotlev1 to be used for making radar chart

# make sure occurrence data for abiotic plot reflects biotic since I cleared abiotic
# trying to say: if occurrence in A is < occurrence for B per plot, then A occurrence = B occurrence
plotlev1$ab_cat <- as.numeric(as.character(plotlev1$ab_cat))

plots <- levels(factor(plotlev1$plot))
grids <- as.numeric(levels(factor(plotlev1$grid)))
spp <- levels(factor(plotlev1$species))

plotlev1$occ <- plotlev1$ab_cat

plotlev1 %>% filter(species == "plaere") -> plotlev2

for(j in grids) {
  plots <- levels(factor(plotlev2$plot[plotlev2$grid == j]))
  for(i in plots) {
    if(length(plotlev2$ab_cat[plotlev2$treatment %in% c("A","B") & plotlev2$plot == i & plotlev2$grid == j]) == 2) {
      if(plotlev2$ab_cat[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j] <= plotlev2$ab_cat[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j]){
        
        plotlev2$occ[plotlev2$treatment == "A" & plotlev2$plot == i & plotlev2$grid == j] <- plotlev2$occ[plotlev2$treatment == "B" & plotlev2$plot == i & plotlev2$grid == j]
      }
    }
  }
}

# run loop again for each species then rbind
plotlev1 %>% filter(species == "miccal") -> plotlev3

for(j in grids) {
  plots <- levels(factor(plotlev3$plot[plotlev3$grid == j]))
  for(i in plots) {
    if(length(plotlev3$ab_cat[plotlev3$treatment %in% c("A","B") & plotlev3$plot == i & plotlev3$grid == j]) == 2) {
      if(plotlev3$ab_cat[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j] <= plotlev3$ab_cat[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j]){
        
        plotlev3$occ[plotlev3$treatment == "A" & plotlev3$plot == i & plotlev3$grid == j] <- plotlev3$occ[plotlev3$treatment == "B" & plotlev3$plot == i & plotlev3$grid == j]
      }
    }
  }
}

plotlev1 %>% filter(species == "brohor") -> plotlev4

for(j in grids) {
  plots <- levels(factor(plotlev4$plot[plotlev4$grid == j]))
  for(i in plots) {
    if(length(plotlev4$ab_cat[plotlev4$treatment %in% c("A","B") & plotlev4$plot == i & plotlev4$grid == j]) == 2) {
      if(plotlev4$ab_cat[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j] <= plotlev4$ab_cat[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j]){
        
        plotlev4$occ[plotlev4$treatment == "A" & plotlev4$plot == i & plotlev4$grid == j] <- plotlev4$occ[plotlev4$treatment == "B" & plotlev4$plot == i & plotlev4$grid == j]
      }
    }
  }
}

plotlev1 %>% filter(species == "vulmic") -> plotlev5

for(j in grids) {
  plots <- levels(factor(plotlev5$plot[plotlev5$grid == j]))
  for(i in plots) {
    if(length(plotlev5$ab_cat[plotlev5$treatment %in% c("A","B") & plotlev5$plot == i & plotlev5$grid == j]) == 2) {
      if(plotlev5$ab_cat[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j] <= plotlev5$ab_cat[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j]){
        
        plotlev5$occ[plotlev5$treatment == "A" & plotlev5$plot == i & plotlev5$grid == j] <- plotlev5$occ[plotlev5$treatment == "B" & plotlev5$plot == i & plotlev5$grid == j]
      }
    }
  }
}

plotlev <- rbind(plotlev2, plotlev3,plotlev4, plotlev5)
rm(plotlev1,plotlev2,plotlev3,plotlev4, plotlev5)
names(plotlev)[11] <- "occurrence" # i know this is spelt wrong but it's ok if i'm consistently wrong ;)

# create y/n
plotlev$occurrence[plotlev$occurrence >= 1] <- c("yes")
plotlev$occurrence[plotlev$occurrence == 0] <- c("no")

plotlev$persistence[plotlev$seed >= 1] <- c("yes")
plotlev$persistence[plotlev$seed == 0] <- c("no")
plotlev$persistence[is.na(plotlev$seed)] <- c("no")

# categories ME, SS_n, SS_y, DL
plotlev$contingency <- with(plotlev, paste0(occurrence, persistence))

plotlev$contingency[plotlev$contingency == "yesyes"] <- c('SS_y')
plotlev$contingency[plotlev$contingency == "nono"] <- c('SS_n')
plotlev$contingency[plotlev$contingency == "yesno"] <- c('ME')
plotlev$contingency[plotlev$contingency == "noyes"] <- c("DL")

# there is a single date in the ab_cat -> fix
plotlev$ab_cat <- as.numeric(plotlev$ab_cat)
# observation 1262 is culprit -> NA = 2
plotlev$ab_cat[is.na(plotlev$ab_cat)] <- 2
#temp_key <- dplyr::select(name_key, -unique, -seed)
plotlev <- left_join(plotlev, nkey, by = c("tag","species","treatment","grid","site"))

# Important note:
# occurrence differing from ab_cat is normal and fine. ab_cat is raw, occurrance considers that abiotic treatment was cleared of naturally occuring species

#### gridlev
# code datasheet for contingencies at grid level
# consider: does taking the average occurrence and persistence measures reflect the population properly??
gridlev1 <- plotlev %>%
  dplyr::select(-11,-12,-13) %>%
  dplyr::group_by(treatment, species, grid) %>%
  dplyr::mutate(grid_ab = mean(ab_cat)) %>%
  dplyr::mutate(grid_n = n()) %>%
  dplyr::mutate(grid_seed = sum(seed)) %>%
  dplyr::select(-tag,-plot,-image,-ab_cat,-status,-seed,-names) %>%
  distinct()
gridlev <- as.data.frame(gridlev1)

gridlev$occurrence <- gridlev$grid_ab
gridlev$persistence[gridlev$grid_seed >= gridlev$grid_n] <- 1
gridlev$persistence[gridlev$grid_seed < gridlev$grid_n] <- 0

gridlev$occurrence[gridlev$occurrence > 0] <- c("yes")
gridlev$occurrence[gridlev$occurrence == 0] <- c("no")
gridlev$persistence[gridlev$persistence >= 1] <- c("yes")
gridlev$persistence[gridlev$persistence < 1] <- c("no")

gridlev$contingency <- with(gridlev, paste0(occurrence, persistence))
gridlev$contingency[gridlev$contingency == "yesyes"] <- c('SS_y')
gridlev$contingency[gridlev$contingency == "nono"] <- c('SS_n')
gridlev$contingency[gridlev$contingency == "yesno"] <- c('ME')
gridlev$contingency[gridlev$contingency == "noyes"] <- c("DL")

rm(gridlev1)

# PROBLEM
# persistence for grid & site = seeds < n = no, seeds >= n = yes

#### sitelev
# contingency yes/no for occurrence and persistence at site level
sitelev1 <- plotlev %>%
  dplyr::select(-11,-12,-13)  %>%
  dplyr::group_by(site, treatment, species) %>%
  dplyr::mutate(site_ab = mean(ab_cat)) %>%
  dplyr::mutate(site_n = n()) %>%
  dplyr::mutate(site_seed = sum(seed)) %>%
  dplyr::select(-tag,-plot,-grid,-image,-ab_cat,-status,-seed,-names) %>%
  distinct()
sitelev <- as.data.frame(sitelev1)

sitelev$occurrence <- sitelev$site_ab
sitelev$persistence[sitelev$site_seed >= sitelev$site_n] <- 1
sitelev$persistence[sitelev$site_seed < sitelev$site_n] <- 0

sitelev$occurrence[sitelev$occurrence > 0] <- c("yes")
sitelev$occurrence[sitelev$occurrence == 0] <- c("no")
sitelev$persistence[sitelev$persistence >= 1] <- c("yes")
sitelev$persistence[sitelev$persistence < 1] <- c("no")

sitelev$contingency <- with(sitelev, paste0(occurrence, persistence))
sitelev$contingency[sitelev$contingency == "yesyes"] <- c('SS_y')
sitelev$contingency[sitelev$contingency == "nono"] <- c('SS_n')
sitelev$contingency[sitelev$contingency == "yesno"] <- c('ME')
sitelev$contingency[sitelev$contingency == "noyes"] <- c("DL")

rm(sitelev1)

#### SAR dt_o dt_p

# Note, in 'plotlev' some ab-cat = 0, occurrence = yes because A was cleared so paired B plot is used to assess occurrence
dt2 <- plotlev %>%
  dplyr::select(tag, species, treatment, grid, site, seed, ab_cat, occurrence, persistence)
dt1 <- left_join(dt2, nkey, by = c("tag", "species", "treatment", "grid", "site"))
dt_p <- dt1 %>%
  dplyr::select(-occurrence, -ab_cat) #%>% dplyr::filter(!persistence%in% 'no')
dt_o <- dt1 %>%
  dplyr::select(-persistence, -seed) #%>% dplyr::filter(!occurrence%in%'no')

dt_p <- left_join(dt_p, dat1, by = "names") # merge with locations data
dt_p <- dplyr::select(dt_p, -tag,-grid,-site)
dt_p$area <- round(dt_p$area, 0)
dt_p <- as.data.frame(dt_p)
dt_p <- distinct(dt_p)

dt_o <- left_join(dt_o, dat1, by = "names") # merge with locations data
dt_o <- dplyr::select(dt_o,-grid,-site)
dt_o$area <- round(dt_o$area, 0)
dt_o <- as.data.frame(dt_o)
dt_o <- distinct(dt_o)

# Bootstrapping data 

## Local
dl <- dplyr::select(plotlev, species, grid, treatment, contingency)

# Abiotic
dla <- dplyr::filter(dl, treatment%in%"A")
dla <- dplyr:: select(dla, -treatment)

# Biotic
# specify treatment
dlb <- dplyr::filter(dl, treatment%in%"B")
dlb <- dplyr:: select(dlb, -treatment)

## Grid

dg <- dplyr::select(gridlev, grid, species, treatment, contingency)

# Abiotic
dga <- dplyr::filter(dg, treatment%in%"A")
dga <- dplyr:: select(dga, -treatment)

# Biotic

dgb <- dplyr::filter(dg, treatment%in%"B")
dgb <- dplyr:: select(dgb, -treatment)

## Site

ds <- dplyr::select(sitelev, site, species, treatment, contingency)

# Abiotic
dsa <- dplyr::filter(ds, treatment%in%"A")
dsa <- dplyr:: select(dsa, -treatment)
dsa <- distinct(dsa)

# Biotic
dsb <- dplyr::filter(ds, treatment%in%"B")
dsb <- dplyr:: select(dsb, -treatment)

# load matrices 
# load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/mat_pb.RData")
# load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/mat_ob.RData")
# load("C:/Users/Megan Szojka/GitHub/SAR-PARs/RData-files/mat_pa.RData")
# checks that all loaded:
# for posting reproducable code start with these files
#View(plotlev)
#View(gridlev) 
#View(sitelev)
#View(df.scales)

#save(plotlev, file = "plotlev.RData")
#save(seeds_clean, file = "seeds_clean.RData")

# MATRICES ####
# ---- don't run, saved! ---------------.
## Persistence
sar_p <- dt_p
sar_p$persistence[sar_p$persistence=='yes'] <- 1
sar_p$persistence[sar_p$persistence == 'no'] <- 0
sar_p$persistence<- as.numeric(sar_p$persistence)

sar_p <- pivot_wider(sar_p, names_from = species, values_from = persistence) # this process converts species columns to a list
sar_p <- dplyr::select(sar_p,-seed)

sar_p$brohor[is.na(sar_p$brohor)] <- 0
sar_p$vulmic[is.na(sar_p$vulmic)] <- 0
sar_p$miccal[is.na(sar_p$miccal)] <- 0
sar_p$plaere[is.na(sar_p$plaere)] <- 0

## Abiotic persistence SAR
# run from here for loop
sar_pa <- sar_p %>%
  dplyr::filter(treatment%in%"A") %>%
  dplyr::select(-treatment)

sar_pa <- dplyr::arrange(sar_pa, run, area)
runs <- c(unique(sar_pa$run))
nrow(sar_pa)
mat_pa <- matrix(NA, ncol = 7, nrow = 801000)
colnames(mat_pa) <- c("names","run","area","plaere","miccal","brohor","vulmic")
k <- 1

for (r in runs){
  temp <- sar_pa %>% dplyr::filter(run%in%r) %>% as.data.frame()
  row <- c(rownames(temp))
  row <- as.numeric(row) # always 1:499
  
  for (i in row){
    ifelse(temp[i-1,4] == 1, temp[i,4]<-1, temp[i,4]<-temp[i,4])
    ifelse(temp[i-1,5] == 1, temp[i,5]<-1, temp[i,5]<-temp[i,5])
    ifelse(temp[i-1,6] == 1, temp[i,6]<-1, temp[i,6]<-temp[i,6])
    ifelse(temp[i-1,7] == 1, temp[i,7]<-1, temp[i,7]<-temp[i,7])
    
    # fill the matrix with current values
    mat_pa[k,1] <- temp[i,1]
    mat_pa[k,2] <- temp[i,2]
    mat_pa[k,3] <- temp[i,3]
    mat_pa[k,4] <- temp[i,4]
    mat_pa[k,5] <- temp[i,5]
    mat_pa[k,6] <- temp[i,6]
    mat_pa[k,7] <- temp[i,7]
    
    k <- k + 1
    
  }
}
mat_pa <- as.data.frame(mat_pa)
class(mat_pa)
class(mat_pa$area)
class(mat_pa$brohor)

mat_pa$brohor<- as.numeric(mat_pa$brohor)
mat_pa$vulmic<- as.numeric(mat_pa$vulmic)
mat_pa$miccal<- as.numeric(mat_pa$miccal)
mat_pa$plaere<- as.numeric(mat_pa$plaere)
mat_pa$area<- as.numeric(mat_pa$area)
mat_pa <- mat_pa %>% 
  mutate(sum =  rowSums(mat_pa[, c(4, 5, 6, 7)]))
class(mat_pa$sum)

## Biotic persistence SAR
# run from here for loop
sar_pb <- sar_p %>%
  dplyr::filter(treatment%in%"B") %>%
  dplyr::select(-treatment)

sar_pb <- dplyr::arrange(sar_pb, run, area)
runs <- c(unique(sar_pb$run))
nrow(sar_pb)
mat_pb <- matrix(NA, ncol = 7, nrow = 727200)
colnames(mat_pb) <- c("names","run","area","plaere","miccal","brohor","vulmic")
k <- 1

for (r in runs){
  temp <- sar_pb %>% dplyr::filter(run%in%r) %>% as.data.frame()
  row <- c(rownames(temp))
  row <- as.numeric(row) # always 1:499
  
  for (i in row){
    ifelse(temp[i-1,4] == 1, temp[i,4]<-1, temp[i,4]<-temp[i,4])
    ifelse(temp[i-1,5] == 1, temp[i,5]<-1, temp[i,5]<-temp[i,5])
    ifelse(temp[i-1,6] == 1, temp[i,6]<-1, temp[i,6]<-temp[i,6])
    ifelse(temp[i-1,7] == 1, temp[i,7]<-1, temp[i,7]<-temp[i,7])
    
    # fill the matrix with current values
    mat_pb[k,1] <- temp[i,1]
    mat_pb[k,2] <- temp[i,2]
    mat_pb[k,3] <- temp[i,3]
    mat_pb[k,4] <- temp[i,4]
    mat_pb[k,5] <- temp[i,5]
    mat_pb[k,6] <- temp[i,6]
    mat_pb[k,7] <- temp[i,7]
    
    k <- k + 1
    
  }
}
mat_pb <- as.data.frame(mat_pb)
class(mat_pb)
class(mat_pb$area)

mat_pb$brohor<- as.numeric(mat_pb$brohor)
mat_pb$vulmic<- as.numeric(mat_pb$vulmic)
mat_pb$miccal<- as.numeric(mat_pb$miccal)
mat_pb$plaere<- as.numeric(mat_pb$plaere)
mat_pb$area<- as.numeric(mat_pb$area)
mat_pb <- mat_pb %>% 
  mutate(sum =  rowSums(mat_pb[, c(4, 5, 6, 7)]))

## Occurrence
sar_o <- dt_o %>%
  dplyr::select(-ab_cat)

sar_o$occurrence[sar_o$occurrence == 'yes'] <- 1
sar_o$occurrence[sar_o$occurrence == 'no'] <- 0
sar_o$occurrence <- as.numeric(sar_o$occurrence)

sar_o <- pivot_wider(sar_o, names_from = species, values_from = occurrence, values_fill = NA)
sar_o$brohor[is.na(sar_o$brohor)] <- 0
sar_o$vulmic[is.na(sar_o$vulmic)] <- 0
sar_o$miccal[is.na(sar_o$miccal)] <- 0
sar_o$plaere[is.na(sar_o$plaere)] <- 0
class(sar_o$brohor)

## Biotic occurrence SAR
# run from here for loop
sar_ob <- sar_o %>%
  dplyr::filter(treatment%in%"B") %>%
  dplyr::select(-treatment,-tag)

sar_ob <- dplyr::arrange(sar_ob, run, area)
runs <- c(unique(sar_ob$run))
nrow(sar_ob) # ONLY LENGTH TO HAVE CHANGED WITH INCLUSION OF HERBIVORY
mat_ob <- matrix(NA, ncol = 7, nrow = 950400)
colnames(mat_ob) <- c("names","run","area","plaere","miccal","brohor","vulmic")
k <- 1

for (r in runs){
  temp <- sar_ob %>% dplyr::filter(run%in%r) %>% as.data.frame()
  row <- c(rownames(temp))
  row <- as.numeric(row) # always 1:499
  
  for (i in row){
    ifelse(temp[i-1,4] == 1, temp[i,4]<-1, temp[i,4]<-temp[i,4])
    ifelse(temp[i-1,5] == 1, temp[i,5]<-1, temp[i,5]<-temp[i,5])
    ifelse(temp[i-1,6] == 1, temp[i,6]<-1, temp[i,6]<-temp[i,6])
    ifelse(temp[i-1,7] == 1, temp[i,7]<-1, temp[i,7]<-temp[i,7])
    
    # fill the matrix with current values
    mat_ob[k,1] <- temp[i,1]
    mat_ob[k,2] <- temp[i,2]
    mat_ob[k,3] <- temp[i,3]
    mat_ob[k,4] <- temp[i,4]
    mat_ob[k,5] <- temp[i,5]
    mat_ob[k,6] <- temp[i,6]
    mat_ob[k,7] <- temp[i,7]
    
    k <- k + 1
    
  }
}
mat_ob <- as.data.frame(mat_ob)
class(mat_ob)
class(mat_ob$area)

mat_ob$brohor<- as.numeric(mat_ob$brohor)
mat_ob$vulmic<- as.numeric(mat_ob$vulmic)
mat_ob$miccal<- as.numeric(mat_ob$miccal)
mat_ob$plaere<- as.numeric(mat_ob$plaere)
mat_ob$area<- as.numeric(mat_ob$area)
mat_ob <- mat_ob %>% 
  mutate(sum =  rowSums(mat_ob[, c(4, 5, 6, 7)]))

# save(mat_pa, file = "mat_pa.RData")
# save(mat_pb, file = "mat_pb.RData")
# save(mat_ob, file = "mat_ob.RData")
#---------------------------------.

# MANUAL FIT ####
# do the 'fit' line as its own model for each facet with (1+log(area)|run) (make sure it's random intercept & slope)

x <-mat_pa%>% dplyr::select(area,run,sum)
mpa <- lm(sum~log(area) + run, data = x) 
plot(fitted(mpa), residuals(mpa))# = fine

coef(mpa)
mpa.c <- coef(mpa)[1]#exp(coef(mpa)[1])
mpa.z <- coef(mpa)[2]
rm(mpa)

x <-mat_pb%>% dplyr::select(area,run,sum)
mpb <- lm(sum~log(area) + run, data = x)
coef(mpb)
mpb.c <- coef(mpb)[1]#exp(coef(mpb)[1])
mpb.z <- coef(mpb)[2]
rm(mpb)

x <-mat_ob%>% dplyr::select(area,run,sum)
mob <- lm(sum ~ log(area) + run, data = x)
coef(mob)
mob.c <- coef(mob)[1]#exp(coef(mob)[1])
mob.z <- coef(mob)[2]
rm(mob)

# calculate y response for lmer models
# v1: using SAR model
# power Model: ln(S) = c + z*ln(A) for log-log 
sar_power = function(x, c, z) {
  c * (x^z) 
}

#logistic model: S = b / (c + A^-z)

#exponential model: S = c + z*ln(A) for semi-log 
sar_exp = function(x, c, z) {
  c + z*log(x)
}

area <- mat_pa %>%
  dplyr::select(area) %>%
  arrange(area) %>%
  distinct()
x <- c(area$area)

resp <- data.frame(area = x, pa = NA, pb = NA, ob = NA)
resp$pa<- sar_exp(x, c=mpa.c,z=mpa.z)
resp$pb<- sar_exp(x, c=mpb.c,z=mpb.z)
resp$ob<- sar_exp(x, c=mob.c,z=mob.z)

#------------------------------------------------------------------------------.
# FIT LINEAR MODEL TO EACH RUN

## mat_pa
runs <- c(unique(mat_pa$run))
coef.pa <- matrix(NA, nrow=900, ncol=3)
colnames(coef.pa) <- c("c","z","run")
k <- 1
for (r in runs){
  s <-mat_pa%>%filter(run == r) %>% dplyr::select(sum)
  x <- mat_pa%>%filter(run == r) %>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area))
  summary(z1)
  coef.pa[k,1] <- coef(z1)[1]
  coef.pa[k,2] <- z1$coefficients[2] 
  coef.pa[k,3] <- r
  k <- k+1
}
coef.pa <- as.data.frame(coef.pa)
coef.pa$c <- as.numeric(coef.pa$c)
coef.pa$z <- as.numeric(coef.pa$z)

## mat_pb
coef.pb <- matrix(NA, nrow=900, ncol=3)
colnames(coef.pb) <- c("c","z","run")
k <- 1
for (r in runs){
  s <-mat_pb%>%filter(run == r)%>%dplyr::select(sum)
  x <- mat_pb%>%filter(run==r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area))
  summary(z1)
  coef.pb[k,1] <- coef(z1)[1]
  coef.pb[k,2] <-  z1$coefficients[2]
  coef.pb[k,3] <- r
  k <- k+1
}
coef.pb <- as.data.frame(coef.pb)
coef.pb$c <- as.numeric(coef.pb$c)
coef.pb$z <- as.numeric(coef.pb$z)

## mat_ob
coef.ob <- matrix(NA, nrow=900, ncol=3)
colnames(coef.ob) <- c("c","z","run")
k <- 1
for (r in runs){
  s<-mat_ob%>%filter(run == r)%>%dplyr::select(sum)
  x<- mat_ob%>%filter(run == r)%>% dplyr::select(area)
  z1 <- lm(s$sum~log(x$area))
  summary(z1)
  coef.ob[k,1] <- coef(z1)[1]
  coef.ob[k,2] <- z1$coefficients[2]
  coef.ob[k,3] <- r
  k <- k+1
}
coef.ob <- as.data.frame(coef.ob)
coef.ob$c <- as.numeric(coef.ob$c)
coef.ob$z <- as.numeric(coef.ob$z)

#------------------------------------------------------------------------------.
# CALCULATE PREDICTED RESPONSES
## responses based on c, A & z

length(x)
length(runs)
y <- data.frame(rep(runs,each=9701), rep(x, times=900))
names(y)[1] <- "run"
names(y)[2] <- "area"
y$resp <- NA
y <- pivot_wider(y,names_from = run, values_from = resp)

# fill in the generic data.frame 'y':
# run in out.pa be column for 'response.pa' and initialize column y before spreading, x is another column repeated for each run
# then call run==r for c and z values in out.pa, fill y data using x 

y_pa <- y
k <- 1
for (i in 2:901){
  y_pa[,i] <- sar_power(x, c=coef.pa[k,1],z=coef.pa[k,2])
  k <- k + 1 # counts the coef.pa row number starting from V1 for each new run
}
y_pb <- y
k <- 1
for (i in 2:901){
  y_pb[,i] <- sar_power(x, c=coef.pb[k,1],z=coef.pb[k,2])
  k <- k+1
}
y_ob <- y
k <- 1
for (i in 2:901){
  y_ob[,i] <- sar_power(x, c=coef.ob[k,1],z=coef.ob[k,2])
  k <- k+1
}

# meld these 3 into one dataframe: sar_all
y_ob1 <- pivot_longer(y_ob, cols = c(2:901), names_to = "run", values_to = "resp")
y_ob1$cat <- "ob"
y_pa1 <- pivot_longer(y_pa, cols = c(2:901), names_to = "run", values_to = "resp")
y_pa1$cat <- "pa"
y_pb1 <- pivot_longer(y_pb, cols = c(2:901), names_to = "run", values_to = "resp")
y_pb1$cat <- "pb"
sar_all <- rbind(y_ob1,y_pa1,y_pb1)
class(sar_all$run)
class(sar_all$cat)
sar_all$run <- as.factor(sar_all$run)
sar_all$cat <- as.factor(sar_all$cat) 


# FIGURE ####
a <- ggplot(data = resp, aes(x=log(area), y=pa)) +
  geom_line()+
  theme_classic() +
  labs(y="Species count", x = "Log area") +
  theme(legend.position = 'none') +
  ylim(0,5)+
  ggtitle("PAR - Fundamental niche") # semi-log model
b <- ggplot(data = resp, aes(x=log(area), y=pb)) +
  geom_line()+
  theme_classic() +
  labs(y="Species count", x = "Log area") +
  theme(legend.position = 'none') +
  ylim(0,5)+
  ggtitle("PAR - Realized niche") # semi-log model
c <- ggplot(data = resp, aes(x=log(area), y=ob)) +
  geom_line()+
  theme_classic() +
  labs(y="Species count", x = "Log area") +
  theme(legend.position = 'none') +
  ylim(0,5)+
  ggtitle("PAR - Realized distribution") # semi-log model
a | b| c
