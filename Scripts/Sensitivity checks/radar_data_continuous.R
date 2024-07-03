##########################################################################.
## This radar plot data uses concept: what how large of scale do we need #
########### to find persistence? Once p = yes, it remains yes ###########
################### at the species levels ###############################
##########################################################################.

# # first change treatment to factor so "A" comes before "B"
# plotlev$treatment<- as.factor(as.character(plotlev$treatment))
# gridlev$treatment<- as.factor(as.character(gridlev$treatment))
# sitelev$treatment<- as.factor(as.character(sitelev$treatment))
# # arrange so that A comes before B
# plotlev <- dplyr::arrange(plotlev,treatment, contingency)
# gridlev <- dplyr::arrange(gridlev,treatment, contingency)
# sitelev <- dplyr::arrange(sitelev,treatment, contingency)

# Use plotlev and if statments to change all 0s to one within a given scale, using matrices_SAR-PAR methods

# plot level use plotlev

# NOTE - not exactly continuous, I still use threshold of grid n < grid seed for persistence
# species level makes it compatable to Figure 3 

# grid level
rad_grid_p <- gridlev
rad_grid_p$persistence[rad_grid_p$persistence=="yes"] <- 1
rad_grid_p$persistence[rad_grid_p$persistence=="no"] <- 0
rad_grid_p$persistence <- as.numeric(rad_grid_p$persistence)

rad_grid_p <- rad_grid_p %>%
  select(treatment, species, grid, persistence) %>%
  group_by(treatment, species, grid) %>%
  pivot_wider(., names_from = "species", values_from = "persistence")
rad_grid_p$brohor[is.na(rad_grid_p$brohor)] <- 0
rad_grid_p$vulmic[is.na(rad_grid_p$vulmic)] <- 0
rad_grid_p$miccal[is.na(rad_grid_p$miccal)] <- 0
rad_grid_p$plaere[is.na(rad_grid_p$plaere)] <- 0
rad_grid_p <- rad_grid_p %>% 
  mutate(sum_persist = rowSums(rad_grid_p[, c(3:6)]))

# repeat for occurrence 
rad_grid_o <- gridlev
rad_grid_o$occurrence[rad_grid_o$occurrence=="yes"] <- 1
rad_grid_o$occurrence[rad_grid_o$occurrence=="no"] <- 0
rad_grid_o$occurrence <- as.numeric(rad_grid_o$occurrence)

rad_grid_o <- rad_grid_o %>%
  select(treatment, species, grid, occurrence) %>%
  group_by(treatment, species, grid) %>%
  pivot_wider(., names_from = "species", values_from = "occurrence")
rad_grid_o$brohor[is.na(rad_grid_o$brohor)] <- 0
rad_grid_o$vulmic[is.na(rad_grid_o$vulmic)] <- 0
rad_grid_o$miccal[is.na(rad_grid_o$miccal)] <- 0
rad_grid_o$plaere[is.na(rad_grid_o$plaere)] <- 0
rad_grid_o <- rad_grid_o %>% 
  mutate(sum_occur = rowSums(rad_grid_o[, c(3:6)]))

# challenge - how to connect specific species and grid in o and p
# pivot longer then join together (unique id issue?)


# site level


# NOW SUMMARIZE FOR THE RADAR PLOTS

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

radat_grid_all <- rad_grid %>%
  ungroup() %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct()
radat_grid_all$ME <- c(0,0)
radat_grid_all$SS_n <- c(0,0)
radat_grid_all = radat_grid_all[,c(1,2,4,5,3)] # rearrrange column order so radar the same

twenty5m2 <- ggradar(radat_grid_all,
                     grid.min = 0,
                     grid.max = 64,
                     legend.position = 'none',
                     background.circle.colour = "white",
                     axis.line.colour = "gray60",
                     gridline.min.colour = "gray60",
                     gridline.mid.colour = "gray60",
                     gridline.max.colour = "gray60",
                     group.colours = cols,
                     grid.label.size = 0)


radat_site_all <- rad_site %>%
  ungroup() %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  dplyr::distinct()
radat_site_all$DL <- c(0,0)
radat_site_all$ME <- c(0,0)
radat_site_all$SS_n <- c(0,0)
radat_site_all = radat_site_all[,c(1,3,4,5,2)] # rearrrange column order so radar the same

one00m2 <- ggradar(radat_site_all,
                   grid.min = 0,
                   grid.max = 24,
                   legend.position = "bottom",
                   background.circle.colour = "white",
                   axis.line.colour = "gray60",
                   gridline.min.colour = "gray60",
                   gridline.mid.colour = "gray60",
                   gridline.max.colour = "gray60",
                   group.colours = cols,
                   grid.label.size = 0) 

library(patchwork)

onem2 +  twenty5m2 + one00m2


#----------------------------------------------------
# DIDN'T WORK
# Plan:
# mat_pa? mat_bp?
# For abiotic comparions, need to join mat_pa and mat_ob by 'name' & 'run' & area' & 'species' in longform
# For biotic comparions, need to join mat_pb and mat_ob by 'name' & 'run' & 'area' & 'species' in longform
# then to both datasets add discrete categories for scales = local, grid, site

# mat 1 and mat2 doing what I want, but occurrence is also continuous version where once = 1 this is true for the rest of the run.
mat_pa1 <- mat_pa %>%
  pivot_longer(., cols = 4:7, names_to = "species", values_to = "persistence_a") %>%
  select(-sum)
mat_ob1 <- mat_ob %>%
  pivot_longer(., cols = 4:7,names_to = "species", values_to = "occurrence") %>%
  select(-sum)
mat_pb1 <- mat_pb %>%
  pivot_longer(., cols = 4:7,names_to = "species", values_to = "persistence_b") %>%
  select(-sum)

mat1 <- full_join(mat_pa1, mat_ob1, by = c("names","run","area", "species"))
mat2 <- full_join(mat_pb1, mat_ob1, by = c("names","run","area", "species"))

mat1[is.na(mat1)] <- 0
mat1$contingency <- with(mat1, paste0(occurrence, persistence_a))
mat1$contingency[mat1$contingency == "11"] <- c('SS_y')
mat1$contingency[mat1$contingency == "00"] <- c('SS_n')
mat1$contingency[mat1$contingency == "10"] <- c('ME')
mat1$contingency[mat1$contingency == "01"] <- c("DL") 

mat2[is.na(mat2)] <- 0
mat2$contingency <- with(mat2, paste0(occurrence, persistence_b))
mat2$contingency[mat2$contingency == "11"] <- c('SS_y')
mat2$contingency[mat2$contingency == "00"] <- c('SS_n')
mat2$contingency[mat2$contingency == "10"] <- c('ME')
mat2$contingency[mat2$contingency == "01"] <- c("DL") 


rm(mat_pa1,mat_ob1,mat_pb1)

## Now summary stats on mat1 and mat2
# still need to add breaks in the areas considered above
# occurrence is not held over within a run, persistence is 
sumary_mat1 <- mat1 %>%
  as.data.frame() %>%
  group_by(run, area, contingency) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::summarize(mean = mean(freq)) %>%#, se = stderr(freq)
  ungroup() %>%
  select(-run) %>%
  distinct()
sumary_mat2 <- mat2 %>%
  as.data.frame() %>%
  group_by(run, area, contingency) %>%
  dplyr::mutate(freq = n()) %>%
  dplyr::summarize(mean = mean(freq)) %>%#, se = stderr(freq)
  ungroup() %>%
  select(-run) %>%
  distinct()

sumary_mat1$discrete_area <- cut(sumary_mat1$area, breaks = 3, labels = c("local","community","regional"))
ggplot(data = sumary_mat1, aes(y = mean, x=discrete_area, fill=contingency)) +
  geom_bar(stat = 'identity') 
sumary_mat2$discrete_area <- cut(sumary_mat2$area, breaks = 3, labels = c("local","community","regional"))
ggplot(data = sumary_mat2, aes(y = mean, x=discrete_area, fill=contingency)) +
  geom_bar(stat = 'identity') 

# Compare with previous way:
#using ...? plotlev?
ggplot(data = filter(plotlev, treatment %in% "A"), aes(y = mean, x=discrete_area, fill=contingency)) +
  geom_bar(stat = 'identity') 


## ALL SPECIES ####

## plot level
radat_plot_all <- plotlev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame() # PROBLEM - with order of A vs B! This dataframe A is row #2

##View(radat_plot_all)

rownames(radat_plot_all) <- c("Abiotic", "Biotic")
# reorder columns
radat_plot_all <- dplyr::select(radat_plot_all,DL,ME,SS_n,SS_y)
colnames(radat_plot_all) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_plot_all)

# use max to denote axis length
total <-rbind(rep(638,4), rep(0,4)) # max=620
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_plot_all <- rbind(total, radat_plot_all)

## gridlev
radat_grid_all <- gridlev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

rownames(radat_grid_all) <- c("Abiotic", "Biotic")

##View(radat_grid_all)

# reorder columns
radat_grid_all <- dplyr::select(radat_grid_all,DL,ME,SS_n,SS_y)
colnames(radat_grid_all) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_grid_all) #44

total <-rbind(rep(44,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_all <- rbind(total, radat_grid_all) #yes

# sitlev
radat_site_all <- sitelev %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::arrange(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  dplyr::distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment)

View(radat_site_all)

# have to create dataframe based on above values, as DL and SS_n are non-existent
DL <- c(0,0)
ME <- c(6,14)
SS_n <- c(0,0)
SS_y <- c(18,10)
radat_site_all <- data.frame (DL,ME,SS_n,SS_y)
colnames(radat_site_all) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
rownames(radat_site_all) <- c("Abiotic","Biotic") 

total <-rbind(rep(18,4), rep(0,4)) #n=14
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_site_all <- rbind(total, radat_site_all)

