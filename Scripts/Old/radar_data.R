##########################################################################.
## Moved all set up for radar plots here - Plotting done in Paper1_figues #
##########################################################################.

# first change treatment to factor so "A" comes before "B"
plotlev$treatment<- as.factor(as.character(plotlev$treatment))
gridlev$treatment<- as.factor(as.character(gridlev$treatment))
sitelev$treatment<- as.factor(as.character(sitelev$treatment))
# arrange so that A comes before B
plotlev <- dplyr::arrange(plotlev,treatment, contingency)
gridlev <- dplyr::arrange(gridlev,treatment, contingency)
sitelev <- dplyr::arrange(sitelev,treatment, contingency)

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

#-----------------------------------------------------------------------.
## BY SPECIES ####

#### miccal ####

## plot level
radat_plot_miccal <- plotlev %>%
  filter(species%in%"miccal") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

##View(radat_plot_miccal)

# reorder columns
radat_plot_miccal <- dplyr::select(radat_plot_miccal,DL,ME,SS_n,SS_y)
rownames(radat_plot_miccal) <- c("Abiotic", "Biotic")
colnames(radat_plot_miccal) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_plot_miccal)
# use max to denote axis length
total <-rbind(rep(215,4), rep(0,4)) # max=211
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_plot_miccal <- rbind(total, radat_plot_miccal)

## gridlev
radat_grid_miccal <- gridlev %>%
  filter(species%in%"miccal") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

View(radat_grid_miccal)

DL <- c(0,0)
ME <- c(15,18)
SS_n <- c(1,0)
SS_y <- c(2,0)
radat_grid_miccal <- data.frame(DL,ME,SS_n,SS_y)
rownames(radat_grid_miccal) <- c("Abiotic", "Biotic")
colnames(radat_grid_miccal) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

total <-rbind(rep(18,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_miccal <- rbind(total, radat_grid_miccal) #yes

# sitlev
radat_site_miccal <- sitelev %>%
  filter(species%in%"miccal") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  dplyr::distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment)

##View(radat_site_miccal)

# have to create dataframe based on above values, as DL and SS_n are non-existent
DL <- c(0,0)
ME <- c(6,6)
SS_n <- c(0,0)
SS_y <- c(0,0)
radat_site_miccal <- data.frame (DL,ME,SS_n,SS_y)
rownames(radat_site_miccal) <- c("Abiotic","Biotic") 
colnames(radat_site_miccal) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

total <-rbind(rep(6,4), rep(0,4))
colnames(total) <-c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_site_miccal <- rbind(total, radat_site_miccal)

#-----------------------------------------------------------------------.
#### plaere ####

## plot level
radat_plot_plaere <- plotlev %>%
  filter(species%in%"plaere") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_plot_plaere)

# reorder columns
radat_plot_plaere <- dplyr::select(radat_plot_plaere,DL,ME,SS_n,SS_y)

rownames(radat_plot_plaere) <- c("Abiotic", "Biotic")
colnames(radat_plot_plaere) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_plot_plaere)
total <-rbind(rep(156,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_plot_plaere <- rbind(total, radat_plot_plaere)

## grid level
radat_grid_plaere <- gridlev %>%
  filter(species%in%"plaere") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_grid_plaere)
# reorder columns
radat_grid_plaere <- dplyr::select(radat_grid_plaere,DL,ME,SS_n,SS_y)

rownames(radat_grid_plaere) <- c("Abiotic", "Biotic")
colnames(radat_grid_plaere) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_plaere[1,1] <- 0
radat_grid_plaere[2,3] <- 0


max(radat_grid_plaere)
total <-rbind(rep(14,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_plaere <- rbind(total, radat_grid_plaere)

## site level
radat_site_plaere <- sitelev %>%
  filter(species%in%"plaere") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_site_plaere)

SS_y <- c(6,3)
ME <- c(0,3)
DL <- c(0,0)
SS_n <- c(0,0)

radat_site_plaere <- data.frame(DL,ME,SS_n,SS_y)
rownames(radat_site_plaere) <- c("Abiotic", "Biotic")
colnames(radat_site_plaere) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_site_plaere)
total <-rbind(rep(6,4), rep(0,4))
colnames(total) <-c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_site_plaere <- rbind(total, radat_site_plaere)

#-----------------------------------------------------------------------.
#### vulmic ####

## plot level
radat_plot_vulmic <- plotlev %>%
  filter(species%in%"vulmic") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_plot_vulmic)

# reorder columns
radat_plot_vulmic <- dplyr::select(radat_plot_vulmic,DL,ME,SS_n,SS_y)

rownames(radat_plot_vulmic) <- c("Abiotic", "Biotic")
colnames(radat_plot_vulmic) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_plot_vulmic)
total <-rbind(rep(143,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_plot_vulmic <- rbind(total, radat_plot_vulmic)


## grid level
radat_grid_vulmic <- gridlev %>%
  filter(species%in%"vulmic") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_grid_vulmic)

radat_grid_vulmic[1,4] <- 0
radat_grid_vulmic[2,2] <- 0

# reorder columns
radat_grid_vulmic <- dplyr::select(radat_grid_vulmic,DL,ME,SS_n,SS_y)

rownames(radat_grid_vulmic) <- c("Abiotic", "Biotic")
colnames(radat_grid_vulmic) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_grid_vulmic)
total <-rbind(rep(16,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_vulmic <- rbind(total, radat_grid_vulmic)


## site level
radat_site_vulmic <- sitelev %>%
  filter(species%in%"vulmic") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_site_vulmic)

SS_y <- c(6,5)
ME <- c(0,1)
DL <- c(0,0)
SS_n <- c(0,0)

radat_site_vulmic <- data.frame(DL, ME, SS_n, SS_y)
rownames(radat_site_vulmic) <- c("Abiotic", "Biotic")
colnames(radat_site_vulmic) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_site_vulmic)
total <-rbind(rep(6,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_site_vulmic <- rbind(total, radat_site_vulmic)

#-----------------------------------------------------------------------.
#### brohor ####

## plot level
radat_plot_brohor <- plotlev %>%
  filter(species%in%"brohor") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_plot_brohor)

# reorder columns
radat_plot_brohor <- dplyr::select(radat_plot_brohor,DL,ME,SS_n,SS_y)

rownames(radat_plot_brohor) <- c("Abiotic", "Biotic")
colnames(radat_plot_brohor) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_plot_brohor)
total <-rbind(rep(165,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_plot_brohor <- rbind(total, radat_plot_brohor)

## grid level
radat_grid_brohor <- gridlev %>%
  filter(species%in%"brohor") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

#View(radat_grid_brohor)

# reorder columns
radat_grid_brohor <- dplyr::select(radat_grid_brohor,DL,ME,SS_n,SS_y)

rownames(radat_grid_brohor) <- c("Abiotic", "Biotic")
colnames(radat_grid_brohor) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_grid_brohor)
total <-rbind(rep(12,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_grid_brohor <- rbind(total, radat_grid_brohor)

## site level
radat_site_brohor <- sitelev %>%
  filter(species%in%"brohor") %>%
  dplyr::select(treatment, contingency) %>%
  dplyr::group_by(treatment, contingency) %>%
  dplyr::add_count() %>%
  distinct() %>%
  pivot_wider(names_from = contingency, values_from = n) %>%
  distinct() %>%
  dplyr::ungroup(treatment) %>%
  dplyr::select(-treatment) %>%
  as.data.frame()

##View(radat_site_brohor)

SS_y <- c(6,2)
ME <- c(0,4)
DL <- c(0,0)
SS_n <- c(0,0)

radat_site_brohor <- data.frame(DL, ME, SS_n, SS_y)
rownames(radat_site_brohor) <- c("Abiotic", "Biotic")
colnames(radat_site_brohor) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")

max(radat_site_brohor)
total <-rbind(rep(6,4), rep(0,4))
colnames(total) <- c("Dispersal limitation","Source-sink","Species sorting (-)","Species sorting (+)")
radat_site_brohor <- rbind(total, radat_site_brohor)

## FIN ##