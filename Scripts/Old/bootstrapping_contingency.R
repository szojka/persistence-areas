################################################################################
#################### Bootstrapping Contingencies ###############################
################################################################################

# Obj 1 ####
# How frequently does realized niches differ from fundamental niches, and what is the role of biotic interactions in changing these?

# Test 1) thinking only of communities (==grids, n=18) …
# use Chi-test of independence (Anova works on means and variance, I have observed n & expected n)

# does the frequency of community statuses differ between biotic and abiotic treatments?

# METHOD 2 ####
# resampling data with bootstraping, telling us what we expect the real distribution to be, and how close it is to null (0.25) 

# use dla & dlb; dga & dgb; dsa & dsb without zeros
# at the end you will have 6 data frames from which I'll make CI round the 'real' data point. Make CI's by length*0.95 = value. Discard min and max observations falling within that value.

# for each run (900)
# for each grid (18)
# for  each species (4)
# sample contingency
# calculate mean and put into final data
nrunz <- 1000 # bump up to 10000 once things work

# Local level + Abiotic

#NaNs bc no brohor or plaere observation in grid 26. 
dla <- filter(dla, !grid %in% 26)

grids <- c(unique(dla$grid))
spp <- c(unique(dla$species))
mat_dla <- matrix(NA, nrow = 17*4*nrunz, ncol = 7) # grids=18-1, species=4, runs=10000
colnames(mat_dla) <- c("run","grid","species", "DL", "ME", "SS_n", "SS_y")
k <- 1

for(r in 1:nrunz){
  for(g in grids){
    for(s in spp){
      temp <- dplyr::filter(dla, grid%in%g, species%in%s)
      samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
      # calculate mean occurrence for each contingency
      temp$frequ <- samp
      temp.dl <- sum(temp$frequ == "DL")/length(samp)
      temp.me <- sum(temp$frequ == "ME")/length(samp)
      temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
      temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
      #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
      
      mat_dla[k,7] <- temp.ssy
      mat_dla[k,6] <- temp.ssn
      mat_dla[k,5] <- temp.me
      mat_dla[k,4] <- temp.dl
      mat_dla[k,3] <- s
      mat_dla[k,2] <- g
      mat_dla[k,1] <- r
      
      k <- k+1
    }
  }
} 

# Local level + Biotic
#NaNs bc no brohor or plaere observation in grid 26. 
dlb <- filter(dlb, !grid %in% 26)

grids <- c(unique(dlb$grid))
spp <- c(unique(dlb$species))
mat_dlb <- matrix(NA, nrow = 17*4*nrunz, ncol = 7) # grids=18, species=4, runs=900
colnames(mat_dlb) <- c("run","grid","species", "DL", "ME", "SS_n", "SS_y")
k <- 1
for(r in 1:nrunz){
  for(g in grids){
    for(s in spp){
      temp <- dplyr::filter(dlb, grid%in%g, species%in%s)
      samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
      # calculate mean occurrence for each contingency
      temp$frequ <- samp
      temp.dl <- sum(temp$frequ == "DL")/length(samp)
      temp.me <- sum(temp$frequ == "ME")/length(samp)
      temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
      temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
      #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
      
      mat_dlb[k,7] <- temp.ssy
      mat_dlb[k,6] <- temp.ssn
      mat_dlb[k,5] <- temp.me
      mat_dlb[k,4] <- temp.dl
      mat_dlb[k,3] <- s
      mat_dlb[k,2] <- g
      mat_dlb[k,1] <- r
      
      k <- k+1
    }
  }
}

# Grid + Abiotic
#dga <- filter(dga, !grid %in% 26)

spp <- c(unique(dga$species))
mat_dga <- matrix(NA, nrow = 4*nrunz, ncol = 6) # species=4, runs=900
colnames(mat_dga) <- c("run","species", "DL", "ME", "SS_n", "SS_y")
k <- 1
for(r in 1:nrunz){
  for(s in spp){
    temp <- dplyr::filter(dga, species%in%s)
    samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
    # calculate mean occurrence for each contingency
    temp$frequ <- samp
    temp.dl <- sum(temp$frequ == "DL")/length(samp)
    temp.me <- sum(temp$frequ == "ME")/length(samp)
    temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
    temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
    #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
    
    mat_dga[k,6] <- temp.ssy
    mat_dga[k,5] <- temp.ssn
    mat_dga[k,4] <- temp.me
    mat_dga[k,3] <- temp.dl
    mat_dga[k,2] <- s
    mat_dga[k,1] <- r
    
    k <- k+1
  }
}

# Grid + Biotic
#dgb <- filter(dgb, !grid %in% 26)

spp <- c(unique(dgb$species))
mat_dgb <- matrix(NA, nrow = 4*nrunz, ncol = 6) # species=4, runs=900
colnames(mat_dgb) <- c("run","species", "DL", "ME", "SS_n", "SS_y")
k <- 1
for(r in 1:nrunz){
  for(s in spp){
    temp <- dplyr::filter(dgb, species%in%s)
    samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
    # calculate mean occurrence for each contingency
    temp$frequ <- samp
    temp.dl <- sum(temp$frequ == "DL")/length(samp)
    temp.me <- sum(temp$frequ == "ME")/length(samp)
    temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
    temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
    #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
    
    mat_dgb[k,6] <- temp.ssy
    mat_dgb[k,5] <- temp.ssn
    mat_dgb[k,4] <- temp.me
    mat_dgb[k,3] <- temp.dl
    mat_dgb[k,2] <- s
    mat_dgb[k,1] <- r
    
    k <- k+1
  }
}

# Site + Abiotic
spp <- c(unique(dsa$species))
mat_dsa <- matrix(NA, nrow = 4*nrunz, ncol = 6) # species=4, runs=900
colnames(mat_dsa) <- c("run","species", "DL", "ME", "SS_n", "SS_y")
k <- 1
for(r in 1:nrunz){
  for(s in spp){
    temp <- dplyr::filter(dsa, species%in%s)
    samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
    # calculate mean occurrence for each contingency
    temp$frequ <- samp
    temp.dl <- sum(temp$frequ == "DL")/length(samp)
    temp.me <- sum(temp$frequ == "ME")/length(samp)
    temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
    temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
    #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
    
    mat_dsa[k,6] <- temp.ssy
    mat_dsa[k,5] <- temp.ssn
    mat_dsa[k,4] <- temp.me
    mat_dsa[k,3] <- temp.dl
    mat_dsa[k,2] <- s
    mat_dsa[k,1] <- r
    
    k <- k+1
  }
}

# Site + Biotic
spp <- c(unique(dsb$species))
mat_dsb <- matrix(NA, nrow = 4*nrunz, ncol = 6) # species=4, runs=900
colnames(mat_dsb) <- c("run","species", "DL", "ME", "SS_n", "SS_y")
k <- 1
for(r in 1:nrunz){
  for(s in spp){
    temp <- dplyr::filter(dsb, species%in%s)
    samp <- c(sample(x = temp$contingency, size = length(temp$contingency), replace = T, prob = NULL)) 
    # calculate mean occurrence for each contingency
    temp$frequ <- samp
    temp.dl <- sum(temp$frequ == "DL")/length(samp)
    temp.me <- sum(temp$frequ == "ME")/length(samp)
    temp.ssn <- sum(temp$frequ == "SS_n")/length(samp)
    temp.ssy <- sum(temp$frequ == "SS_y")/length(samp)
    #sum(temp.dl,temp.me,temp.ssn,temp.ssy) = 1 good!
    
    mat_dsb[k,6] <- temp.ssy
    mat_dsb[k,5] <- temp.ssn
    mat_dsb[k,4] <- temp.me
    mat_dsb[k,3] <- temp.dl
    mat_dsb[k,2] <- s
    mat_dsb[k,1] <- r
    
    k <- k+1
  }
}

#----------------------------------------------------------------------------la
# I have data frame, now plot against hline = 0.25
mat_dla <- data.frame(mat_dla)
mat_dla$DL <- as.numeric(mat_dla$DL)
mat_dla$ME <- as.numeric(mat_dla$ME)
mat_dla$SS_n <- as.numeric(mat_dla$SS_n)
mat_dla$SS_y <- as.numeric(mat_dla$SS_y)
mat_dla <- na.omit(mat_dla)

# remove lowest 2.5% of values and highest 2.5% set of values
all_dl <- dplyr::select(mat_dla,DL) 
all_dl <- na.omit(all_dl) 
all_dl <- all_dl %>% arrange(DL)
remove <- length(all_dl$DL)*0.025
all_dl <- head(all_dl, -remove) # this # is length(all_dl$DL)*0.025
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dla,ME) 
all_me <- na.omit(all_me)
remove <- length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dla,SS_n) 
all_ssn <- na.omit(all_ssn)
all_ssn <- all_ssn %>% arrange(SS_n)
remove <- length(all_ssn$SS_n)*0.025
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dla,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_la <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dla$DL),mean(mat_dla$ME),mean(mat_dla$SS_n),mean(mat_dla$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dla$contingency == "DL")/length(dla$contingency),
                                 sum(dla$contingency == "ME")/length(dla$contingency),
                                 sum(dla$contingency == "SS_n")/length(dla$contingency),
                                 sum(dla$contingency == "SS_y")/length(dla$contingency)))

save(avg_la,file = "avg_la.RData") # move loops and avg_X data frames to their own script after

# resample distribution 10,000 times? #### 
# not yet
#----------------------------------------------------------------------------lb
mat_dlb <- data.frame(mat_dlb)
mat_dlb$DL <- as.numeric(mat_dlb$DL)
mat_dlb$ME <- as.numeric(mat_dlb$ME)
mat_dlb$SS_n <- as.numeric(mat_dlb$SS_n)
mat_dlb$SS_y <- as.numeric(mat_dlb$SS_y)
mat_dlb <- na.omit(mat_dlb)

# remove lowest 2.5% of values and highest 2.5% set of values
all_dl <- dplyr::select(mat_dlb,DL) 
all_dl <- na.omit(all_dl)
remove <- length(all_dl$DL)*0.025
all_dl <- all_dl %>% arrange(DL)
all_dl <- head(all_dl, -remove)
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dlb,ME) 
all_me <- na.omit(all_me)
remove <-length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dlb,SS_n) 
all_ssn <- na.omit(all_ssn)
remove <-length(all_ssn$SS_n)*0.025
all_ssn <- all_ssn %>% arrange(SS_n)
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dlb,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_lb <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dlb$DL),mean(mat_dlb$ME),mean(mat_dlb$SS_n),mean(mat_dlb$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dlb$contingency == "DL")/length(dlb$contingency),
                                 sum(dlb$contingency == "ME")/length(dlb$contingency),
                                 sum(dlb$contingency == "SS_n")/length(dlb$contingency),
                                 sum(dlb$contingency == "SS_y")/length(dlb$contingency))) 

save(avg_lb, file = "avg_lb.RData")
#-------------------------------------------------------------------------------ga

#mat_dga 
mat_dga <- data.frame(mat_dga)
mat_dga$DL <- as.numeric(mat_dga$DL)
mat_dga$ME <- as.numeric(mat_dga$ME)
mat_dga$SS_n <- as.numeric(mat_dga$SS_n)
mat_dga$SS_y <- as.numeric(mat_dga$SS_y)
mat_dga <- na.omit(mat_dga)

# remove lowest 2.5% of values and highest 2.5% set of values
all_dl <- dplyr::select(mat_dga,DL) 
all_dl <- na.omit(all_dl)
remove <- length(all_dl$DL)*0.025
all_dl <- all_dl %>% arrange(DL)
all_dl <- head(all_dl, -remove)
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dga,ME) 
all_me <- na.omit(all_me)
remove <- length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dga,SS_n) 
all_ssn <- na.omit(all_ssn)
remove <- length(all_ssn$SS_n)*0.025
all_ssn <- all_ssn %>% arrange(SS_n)
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dga,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_ga <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dga$DL),mean(mat_dga$ME),mean(mat_dga$SS_n),mean(mat_dga$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dga$contingency == "DL")/length(dga$contingency),
                                 sum(dga$contingency == "ME")/length(dga$contingency),
                                 sum(dga$contingency == "SS_n")/length(dga$contingency),
                                 sum(dga$contingency == "SS_y")/length(dga$contingency)))

save(avg_ga,file = "avg_ga.RData") 
#-------------------------------------------------------------------------------gb

#mat_dgb 
mat_dgb <- data.frame(mat_dgb)
mat_dgb$DL <- as.numeric(mat_dgb$DL)
mat_dgb$ME <- as.numeric(mat_dgb$ME)
mat_dgb$SS_n <- as.numeric(mat_dgb$SS_n)
mat_dgb$SS_y <- as.numeric(mat_dgb$SS_y)
mat_dgb <- na.omit(mat_dgb)

all_dl <- dplyr::select(mat_dgb,DL) 
all_dl <- na.omit(all_dl)
remove <- length(all_dl$DL)*0.025
all_dl <- all_dl %>% arrange(DL)
all_dl <- head(all_dl, -remove)
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dgb,ME) 
all_me <- na.omit(all_me)
remove <- length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dgb,SS_n) 
all_ssn <- na.omit(all_ssn)
remove <- length(all_ssn$SS_n)*0.025
all_ssn <- all_ssn %>% arrange(SS_n)
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dgb,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_gb <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dgb$DL),mean(mat_dgb$ME),mean(mat_dgb$SS_n),mean(mat_dgb$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dgb$contingency == "DL")/length(dgb$contingency),
                                 sum(dgb$contingency == "ME")/length(dgb$contingency),
                                 sum(dgb$contingency == "SS_n")/length(dgb$contingency),
                                 sum(dgb$contingency == "SS_y")/length(dgb$contingency)))

save(avg_gb,file = "avg_gb.RData") 
#-------------------------------------------------------------------------------sa

#mat_dsa 
mat_dsa <- data.frame(mat_dsa)
mat_dsa$DL <- as.numeric(mat_dsa$DL)
mat_dsa$ME <- as.numeric(mat_dsa$ME)
mat_dsa$SS_n <- as.numeric(mat_dsa$SS_n)
mat_dsa$SS_y <- as.numeric(mat_dsa$SS_y)
mat_dsa <- na.omit(mat_dsa)

all_dl <- dplyr::select(mat_dsa,DL) 
all_dl <- na.omit(all_dl)
remove <- length(all_dl$DL)*0.025
all_dl <- all_dl %>% arrange(DL)
all_dl <- head(all_dl, -remove)
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dsa,ME) 
all_me <- na.omit(all_me)
remove <- length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dsa,SS_n) 
all_ssn <- na.omit(all_ssn)
remove <- length(all_ssn$SS_n)*0.025
all_ssn <- all_ssn %>% arrange(SS_n)
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dsa,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_sa <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dsa$DL),mean(mat_dsa$ME),mean(mat_dsa$SS_n),mean(mat_dsa$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dsa$contingency == "DL")/length(dsa$contingency),
                                 sum(dsa$contingency == "ME")/length(dsa$contingency),
                                 sum(dsa$contingency == "SS_n")/length(dsa$contingency),
                                 sum(dsa$contingency == "SS_y")/length(dsa$contingency)))

save(avg_sa,file = "avg_sa.RData") 
#-------------------------------------------------------------------------------sb
#mat_dsb 
mat_dsb <- data.frame(mat_dsb)
mat_dsb$DL <- as.numeric(mat_dsb$DL)
mat_dsb$ME <- as.numeric(mat_dsb$ME)
mat_dsb$SS_n <- as.numeric(mat_dsb$SS_n)
mat_dsb$SS_y <- as.numeric(mat_dsb$SS_y)
mat_dsb <- na.omit(mat_dsb)

all_dl <- dplyr::select(mat_dsb,DL) 
all_dl <- na.omit(all_dl)
remove <- length(all_dl$DL)*0.025
all_dl <- all_dl %>% arrange(DL)
all_dl <- head(all_dl, -remove)
all_dl <- tail(all_dl, -remove)
all_dl <- as.data.frame(all_dl)

all_me <- dplyr::select(mat_dsb,ME) 
all_me <- na.omit(all_me)
remove <- length(all_me$ME)*0.025
all_me <- all_me %>% arrange(ME)
all_me <- head(all_me, -remove)
all_me <- tail(all_me, -remove)
all_me <- as.data.frame(all_me)

all_ssn <- dplyr::select(mat_dsb,SS_n) 
all_ssn <- na.omit(all_ssn)
remove <- length(all_ssn$SS_n)*0.025
all_ssn <- all_ssn %>% arrange(SS_n)
all_ssn <- head(all_ssn, -remove)
all_ssn <- tail(all_ssn, -remove)
all_ssn <- as.data.frame(all_ssn)

all_ssy <- dplyr::select(mat_dsb,SS_y) 
all_ssy <- na.omit(all_ssy)
remove <- length(all_ssy$SS_y)*0.025
all_ssy <- all_ssy %>% arrange(SS_y)
all_ssy <- head(all_ssy, -remove)
all_ssy <- tail(all_ssy, -remove)
all_ssy <- as.data.frame(all_ssy)

##Means
avg_sb <- data.frame(contingency=c('DL','ME','SS_n','SS_y'),
                     samp.prop=c(mean(mat_dsb$DL),mean(mat_dsb$ME),mean(mat_dsb$SS_n),mean(mat_dsb$SS_y)),
                     low = c(min(all_dl$DL),min(all_me$ME),min(all_ssn$SS_n),min(all_ssy$SS_y)),
                     high = c(max(all_dl$DL),max(all_me$ME),max(all_ssn$SS_n),max(all_ssy$SS_y)),
                     true.prop=c(sum(dsb$contingency == "DL")/length(dsb$contingency),
                                 sum(dsb$contingency == "ME")/length(dsb$contingency),
                                 sum(dsb$contingency == "SS_n")/length(dsb$contingency),
                                 sum(dsb$contingency == "SS_y")/length(dsb$contingency)))
save(avg_sb,file = "avg_sb.RData") 

# put into one dataframe for summary stats
avg_la$name <- c('la')
avg_lb$name <- c('lb')
avg_ga$name <- c('ga')
avg_gb$name <- c('gb')
avg_sa$name <- c('sa')
avg_sb$name <- c('sb')
avg_all <- rbind(avg_la,avg_lb,avg_ga,avg_gb,avg_sa,avg_sb)
save(avg_all,file = "avg_all.RData") 
write_csv(avg_all, file = "boostrapping summary stats.csv")
## BACKUP Dirichlet Reg 
# http://r-statistics.co/Dirichlet-Regression-With-R.html