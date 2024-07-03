

# ---- don't run, saved! ----------

#-----------------------------------------------------------------------------.
# DESCRIPTION: Wrangle matrices to calculate Sar and Par accumulation off of
#-----------------------------------------------------------------------------.

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
source("Scripts/Source - MAIN fitnessdata.R")

# using dt_o or dt_p

# for each run, organize data by area, when a 1 appears in the species column, it now has a one for each area until the end of the run

## Peristence AND Occupancy (added 3/7/2022 after Rachel lab meeting)
dt_op <- left_join(dt_o, dt_p, by = c('species','treatment','names','run','area'))
dt_op <- dt_op %>% dplyr::filter(occurrence == 'yes') %>% dplyr::select(-ab_cat, -occurrence,-tag)
sar_po <- dt_op
sar_po$persistence[sar_po$persistence=='yes'] <- 1
sar_po$persistence[sar_po$persistence == 'no'] <- 0
sar_po$persistence<- as.numeric(sar_po$persistence)

sar_po <- sar_po %>% # solves unique identifyer problem with pivots wider
  dplyr::group_by(treatment, species) %>%
  dplyr::mutate(row = row_number()) %>%
  tidyr::pivot_wider(names_from = species, values_from = persistence) %>%
  dplyr::select(-row) %>%
  dplyr::ungroup()

sar_po <- dplyr::select(sar_po,-seed)
sar_po
sar_po$brohor[is.na(sar_po$brohor)] <- 0
sar_po$vulmic[is.na(sar_po$vulmic)] <- 0
sar_po$miccal[is.na(sar_po$miccal)] <- 0
sar_po$plaere[is.na(sar_po$plaere)] <- 0

# Abiotic occurrence AND persistence

sar_poa <- sar_po %>%
  dplyr::filter(treatment%in%"A") %>%
  dplyr::select(-treatment)

sar_poa <- dplyr::arrange(sar_poa, run, area)
runs <- c(unique(sar_poa$run))
nrow_sar_poa <- nrow(sar_poa)
mat_poa <- matrix(NA, ncol = 7, nrow = nrow_sar_poa)
colnames(mat_poa) <- c("names","run","area","plaere","miccal","brohor","vulmic")
k <- 1

for (r in runs){
  temp <- sar_poa %>% dplyr::filter(run%in%r) %>% as.data.frame()
  row <- c(rownames(temp))
  row <- as.numeric(row) # always 1:499
  
  for (i in row){
    ifelse(temp[i-1,4] == 1, temp[i,4]<-1, temp[i,4]<-temp[i,4])
    ifelse(temp[i-1,5] == 1, temp[i,5]<-1, temp[i,5]<-temp[i,5])
    ifelse(temp[i-1,6] == 1, temp[i,6]<-1, temp[i,6]<-temp[i,6])
    ifelse(temp[i-1,7] == 1, temp[i,7]<-1, temp[i,7]<-temp[i,7])
    
    # fill the matrix with current values
    mat_poa[k,1] <- temp[i,1]
    mat_poa[k,2] <- temp[i,2]
    mat_poa[k,3] <- temp[i,3]
    mat_poa[k,4] <- temp[i,4]
    mat_poa[k,5] <- temp[i,5]
    mat_poa[k,6] <- temp[i,6]
    mat_poa[k,7] <- temp[i,7]
    
    k <- k + 1
    
  }
}
mat_poa <- as.data.frame(mat_poa)
class(mat_poa)
class(mat_poa$area)

mat_poa$brohor<- as.numeric(mat_poa$brohor)
mat_poa$vulmic<- as.numeric(mat_poa$vulmic)
mat_poa$miccal<- as.numeric(mat_poa$miccal)
mat_poa$plaere<- as.numeric(mat_poa$plaere)
mat_poa$area<- as.numeric(mat_poa$area)
mat_poa <- mat_poa %>% 
  mutate(sum =  rowSums(mat_poa[, c(4, 5, 6, 7)]))

# Biotic occurrence AND persistence 
sar_pob <- sar_po %>%
  dplyr::filter(treatment%in%"B") %>%
  dplyr::select(-treatment)

sar_pob <- dplyr::arrange(sar_pob, run, area)
runs <- c(unique(sar_pob$run))
nrow_sar_pob <- nrow(sar_pob)
mat_pob <- matrix(NA, ncol = 7, nrow = nrow_sar_pob)
colnames(mat_pob) <- c("names","run","area","plaere","miccal","brohor","vulmic")
k <- 1

for (r in runs){
  temp <- sar_pob %>% dplyr::filter(run%in%r) %>% as.data.frame()
  row <- c(rownames(temp))
  row <- as.numeric(row) # always 1:499
  
  for (i in row){
    ifelse(temp[i-1,4] == 1, temp[i,4]<-1, temp[i,4]<-temp[i,4])
    ifelse(temp[i-1,5] == 1, temp[i,5]<-1, temp[i,5]<-temp[i,5])
    ifelse(temp[i-1,6] == 1, temp[i,6]<-1, temp[i,6]<-temp[i,6])
    ifelse(temp[i-1,7] == 1, temp[i,7]<-1, temp[i,7]<-temp[i,7])
    
    # fill the matrix with current values
    mat_pob[k,1] <- temp[i,1]
    mat_pob[k,2] <- temp[i,2]
    mat_pob[k,3] <- temp[i,3]
    mat_pob[k,4] <- temp[i,4]
    mat_pob[k,5] <- temp[i,5]
    mat_pob[k,6] <- temp[i,6]
    mat_pob[k,7] <- temp[i,7]
    
    k <- k + 1
    
  }
}
mat_pob <- as.data.frame(mat_pob)
class(mat_pob)
class(mat_pob$area)

mat_pob$brohor<- as.numeric(mat_pob$brohor)
mat_pob$vulmic<- as.numeric(mat_pob$vulmic)
mat_pob$miccal<- as.numeric(mat_pob$miccal)
mat_pob$plaere<- as.numeric(mat_pob$plaere)
mat_pob$area<- as.numeric(mat_pob$area)
mat_pob <- mat_pob %>% 
  mutate(sum =  rowSums(mat_pob[, c(4, 5, 6, 7)]))


## Persistence
sar_p <- dt_p
sar_p$persistence[sar_p$persistence=='yes'] <- 1
sar_p$persistence[sar_p$persistence == 'no'] <- 0
sar_p$persistence<- as.numeric(sar_p$persistence)

sar_p <- pivot_wider(sar_p, names_from = species, values_from = persistence) 
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
nrow_sar_pa <- nrow(sar_pa)
mat_pa <- matrix(NA, ncol = 7, nrow = nrow_sar_pa)
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
nrow_sar_pb <- nrow(sar_pb)
mat_pb <- matrix(NA, ncol = 7, nrow = nrow_sar_pb)
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
  dplyr::select(-treatment, -tag)

sar_ob <- dplyr::arrange(sar_ob, run, area)
runs <- c(unique(sar_ob$run))
nrow_sar_ob <- nrow(sar_ob)
mat_ob <- matrix(NA, ncol = 7, nrow = nrow_sar_ob)
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

# eleminate repeats in matrices (shouldn't be more than one obs by run, area, name)
mat_pa <- distinct(mat_pa) 
mat_pb <- distinct(mat_pb) 
mat_ob <- distinct(mat_ob) 
mat_poa <- distinct(mat_poa)
mat_pob <- distinct(mat_pob)


save(mat_pa, file = "mat_pa.RData")
save(mat_pb, file = "mat_pb.RData")
save(mat_ob, file = "mat_ob.RData")
save(mat_poa, file = "mat_poa.RData")
save(mat_pob, file = "mat_pob.RData")
#---------------------------------.