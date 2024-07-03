#  Much of this is still in Source_fitness data, but 'perf' data etc had been removed. Find it below it needed.

dl <- dplyr::select(plotlev, species, grid, treatment, contingency)

# perfect data
species <- c('brohor', 'miccal', 'plaere', 'vulmic')
contingency <- c("DL","ME","SS_n","SS_y")
length(unique(names))
perf <- data.frame(names = as.factor(rep(names,each=4)),
                   species = as.factor(rep(species,times=4)),
                   contingency = as.factor(rep(contingency, each=1800)))
perf$code <- paste(perf$names, perf$species, sep="_")
perf$unique <- paste(perf$names,perf$species,perf$contingency, sep="_")

## Abiotic
dla <- dplyr::filter(dl, treatment%in%"A")
dla <- dplyr:: select(dla, -treatment)
dla$data <- 1

dla_new <- left_join(perf, dla, by = c("names","species","contingency"))

dla_new$data<-ifelse(is.na(dla_new$data)==TRUE,0,1) # contingencies aren't working

Sum <- c()
k <- 1
for (i in unique(dla_new$code)){
  temp <- dplyr::filter(dla_new, code == i)
  Sum[k] <- sum(temp$data)
  k <- k+1
}
Code <- unique(dla_new$code)
filt <- data.frame(Code, Sum)
filt <- filter(filt, Sum==0)

perf_la <- dplyr::filter(perf, !code%in%filt$Code)
nsg <- dla %>%
  dplyr::select(names, grid, site) %>%
  distinct()
perf_la <- left_join(nsg, perf_la, by = "names")

dla <- left_join(perf_la, dla, by = c("names","species","contingency","grid","site"))
dla$grid <- as.factor(as.character(dla$grid))

dla$data<-ifelse(is.na(dla$data)==TRUE,0,1) # contingencies aren't working
#check <- dplyr::filter(dla, names%in%"1_HUT11")

## Biotic
# specify treatment
dlb <- dplyr::filter(dl, treatment%in%"B")
dlb <- dplyr:: select(dlb, -treatment)
dlb$data <- 1

dlb_new <- left_join(perf, dlb, by = c("names","species","contingency"))

dlb_new$data<-ifelse(is.na(dlb_new$data)==TRUE,0,1) # contingencies aren't working

# take that list of removed codes and alter perf data accordingly
Sum <- c()
k <- 1
for (i in unique(dlb_new$code)){
  temp <- dplyr::filter(dlb_new, code == i)
  Sum[k] <- sum(temp$data)
  k <- k+1
}
Code <- unique(dlb_new$code)
filt <- data.frame(Code, Sum)
filt <- filter(filt, Sum==0)

perf_lb <- dplyr::filter(perf, !code%in%filt$Code)
nsg <- dlb %>%
  dplyr::select(names, grid, site) %>%
  distinct()
perf_lb <- left_join(nsg, perf_lb, by = "names")

dlb <- left_join(perf_lb, dlb, by = c("names","species","contingency","grid","site"))
dlb$grid <- as.factor(as.character(dlb$grid))

dlb$data<-ifelse(is.na(dlb$data)==TRUE,0,1) # contingencies aren't working
check <- dplyr::filter(dlb, names%in%"1_HUT11")

#### regional level ####
site_names <- plotlev %>%
  dplyr::select(site, grid) %>%
  distinct()
class(site_names$grid)
site_names$grid <- as.factor(as.character(site_names$grid))

dg <- dplyr::select(gridlev, grid, species, treatment, contingency)

## Abiotic
dga <- dplyr::filter(dg, treatment%in%"A")
dga <- dplyr:: select(dga, -treatment)
dga$data <- 1

# perfect data
species <- c('brohor', 'miccal', 'plaere', 'vulmic')
contingency <- c("DL","ME","SS_n","SS_y")
grid_names <- unique(dga$grid)
perf <- data.frame(grid = as.factor(rep(grid_names,each=4)),
                   species = as.factor(rep(species,times=4)),
                   contingency = as.factor(rep(contingency, each=72)))
perf$code <- paste(perf$grid, perf$species, sep="_")
perf$unique <- paste(perf$grid,perf$species,perf$contingency, sep="_")
dga$grid <- as.factor(as.character(dga$grid))

dga_new <- left_join(perf, dga, by = c("grid","species","contingency"))

dga_new$data<-ifelse(is.na(dga_new$data)==TRUE,0,1)

# take that list of removed codes and alter perf data accordingly
Sum <- c()
k <- 1
for (i in unique(dga_new$code)){
  temp <- dplyr::filter(dga_new, code == i)
  Sum[k] <- sum(temp$data)
  k <- k+1
}
Code <- unique(dga_new$code)
filt <- data.frame(Code, Sum)
filt <- filter(filt, Sum==0)

perf_ga <- dplyr::filter(perf, !code%in%filt$Code)
dga <- left_join(perf_ga, dga, by = c("grid","species","contingency"))
dga$data<-ifelse(is.na(dga$data)==TRUE,0,1)
dga <- distinct(dga)
check <- dplyr::filter(dga, grid%in%"1")
dga <- left_join(dga, site_names, by = "grid")

#### site level ####

ds <- dplyr::select(sitelev, site, species, treatment, contingency)

# specify treatment
dsa <- dplyr::filter(ds, treatment%in%"A")
dsa <- dplyr:: select(dsa, -treatment)
# dsa$data <- 1
dsa <- distinct(dsa)

# perf data for each scale will look diff, need unique identifier
# perfect data
species <- c('brohor', 'miccal', 'plaere', 'vulmic')
contingency <- c("DL","ME","SS_n","SS_y")
sites <- unique(dsa$site)
perf <- data.frame(site = as.factor(rep(sites,each=4)),
                   species = as.factor(rep(species,times=4)),
                   contingency = as.factor(rep(contingency, each=24)))
perf$code <- paste(perf$site, perf$species, sep="_")
perf$unique <- paste(perf$site,perf$species,perf$contingency, sep="_")

class(perf$site)
class(dsa$site)
dsa$site <- as.factor(as.character(dsa$site))

dsa_new <- left_join(perf, dsa, by = c("site","species","contingency"))

dsa_new$data<-ifelse(is.na(dsa_new$data)==TRUE,0,1)

# take that list of removed codes and alter perf data accordingly
Sum <- c()
k <- 1
for (i in unique(dsa_new$code)){
  temp <- dplyr::filter(dsa_new, code == i)
  Sum[k] <- sum(temp$data)
  k <- k+1
}
Code <- unique(dsa_new$code)
filt <- data.frame(Code, Sum)
filt <- filter(filt, Sum==0)

perf_sa <- dplyr::filter(perf, !code%in%filt$Code)
dsa <- left_join(perf_sa, dsa, by = c("site","species","contingency"))
dsa$data<-ifelse(is.na(dsa$data)==TRUE,0,1)
check <- dplyr::filter(dsa, site%in%"TAIL")

# Missing dgb & dsb
