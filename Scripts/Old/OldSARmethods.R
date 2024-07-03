# Old SAR fitting methods
~~~~~~~~~~~~~~~~~~~~~~~~~


## METHOD 1: specaccum & fitspecaccum from vegan package

# Storch 2016: "By ‘sample areas’ I mean contiguous plots of given area – it makes no sense to speak about the SAR if the areas are not contiguous (Dengler 2008, 2009; Dengler & Oldeland 2010; Güler et al. 2016), i.e. in the case of so‐called species accumulation curves or collector curves. Such curves are driven by other factors (namely pure sampling effects) and their behaviour is thus different"
# therefore not a good method for SAR :(

# THE NEXT TWO METHODS DEAL WITH SPECIES AREA RELATIONSHIP SPECIFICALLY. IN CHOOSING WHICH TO USE CONSIDER THAT MY DATA IS NESTED, AND NOT LIKE ISLANDS
# S = cA^z or S=kA^z
# for methods see https://esajournals.onlinelibrary.wiley.com/doi/10.1002/ecy.3231

## METHOD 2: SSarrhenius from vegan package (self-start models)
# because species accumulation curves don't = species-area relationship...
# functions: SSarrhenius / SSgleason/ SSgitay/ SSlomolino
##### PROBLEMS
# works but plot is v odd, because at large scales 1 still shows up..
# not summing species as area increases (difference bw SAR and SAC?)
# problem also is that data is arranged as islands, not nested

## METHOD 3: mmSAR package
#install.packages("mmSAR", repos="http://R-Forge.R-project.org")
# Example:
#loading all available models
# data(power) # arrhenisu
# data(expo) #gleason, convex, not asymptotic
# data(negexpo)
# data(monod)
# data(ratio)
# data(logist)
# data(lomolino)
# data(weibull)
# #loading the Galapagos Islands plants data set (Preston, 1962)
# data(data.arr)
# data("data.galap")
# #creating a vector of model names
# mods = c("power","expo","negexpo","monod","logist","ratio","lomolino","weibull")
# #fitting all the models to the Galapagos dataset and perform multimodel averaging
# resAverage <- multiSAR(modelList=mods,data.arr)

## METHOD 4 (final): raw data ####
# using matrices mat_oa, mat_ob etc.

# ---- don't run, saved! ----------
library(mmSAR)

data(power) # might be creating problems since I have finite data
data(negexpo) # convex asymptotic: s == c * (1 - exp(-z * a))
?negexpo

mmPA1 <- mat_pa %>%
  dplyr::select(run,area,sum)
# make a loop to get these values for all runs, save values in matrix.
out.pa <- matrix(data = NA, nrow = 900, ncol = 4)
colnames(out.pa) <- c("run",'c','z','r2')
runs <- unique(mmPA1$run)
k <- 1
for (r in runs){
  temp <- subset(mmPA1,run == r)
  temp <- dplyr::select(temp,-run) # creat mmSAR object
  names(temp)[1] <- 'a'
  names(temp)[2] <- 's'
  name = 'Arrhenius 1921'
  data = as.data.frame(temp)
  mmPA <- list(name, data)
  names(mmPA) <- c("name","data")
  
  output <- rssoptim(negexpo, mmPA, norTest = "lillie", verb=F) # fit run to SAR negexpo model
  #print(output$par[1])
  
  # save values to dataframe
  out.pa[k,1] <- r
  out.pa[k,2] <- output$par[1]
  out.pa[k,3] <- output$par[2]
  out.pa[k,4] <- output$R2
  
  k <- k+1
}
View(out.pa)
out.pa <- as.data.frame(out.pa)
out.pa$c <- as.numeric(out.pa$c)
out.pa$z <- as.numeric(out.pa$z)

mmPB1 <- mat_pb %>%
  dplyr::select(run,area,sum)
# make a loop to get these values for all runs, save values in matrix.
out.pb <- matrix(data = NA, nrow = 900, ncol = 4)
colnames(out.pb) <- c("run",'c','z','r2')
runs <- unique(mmPB1$run)
k <- 1
for (r in runs){
  temp <- subset(mmPB1,run == r)
  temp <- dplyr::select(temp,-run) # creat mmSAR object
  names(temp)[1] <- 'a'
  names(temp)[2] <- 's'
  name = 'Arrhenius 1921'
  data = as.data.frame(temp)
  mmPB <- list(name, data)
  names(mmPB) <- c("name","data")
  
  output <- rssoptim(negexpo, mmPB, norTest = "lillie", verb=F) # fit run to SAR negexpo model
  #print(output$par[1])
  
  # save values to dataframe
  out.pb[k,1] <- r
  out.pb[k,2] <- output$par[1]
  out.pb[k,3] <- output$par[2]
  out.pb[k,4] <- output$R2
  
  k <- k+1
}
View(out.pb)
out.pb <- as.data.frame(out.pb)
out.pb$c <- as.numeric(out.pb$c)
out.pb$z <- as.numeric(out.pb$z)

mmOB1 <- mat_ob %>%
  dplyr::select(run,area,sum)
# make a loop to get these values for all runs, save values in matrix.
out.ob <- matrix(data = NA, nrow = 900, ncol = 4)
colnames(out.ob) <- c("run",'c','z','r2')
runs <- unique(mmOB1$run)
k <- 1
for (r in runs){
  temp <- subset(mmOB1,run == r)
  temp <- dplyr::select(temp,-run) # create mmSAR object
  names(temp)[1] <- 'a'
  names(temp)[2] <- 's'
  name = 'Arrhenius 1921'
  data = as.data.frame(temp)
  mmPB <- list(name, data)
  names(mmPB) <- c("name","data")
  
  output <- rssoptim(negexpo, mmPB, norTest = "lillie", verb=F) # fit run to SAR negexpo model
  #print(output$par[1])
  
  # save values to dataframe
  out.ob[k,1] <- r
  out.ob[k,2] <- output$par[1]
  out.ob[k,3] <- output$par[2]
  out.ob[k,4] <- output$R2
  
  k <- k+1
}
View(out.ob)
out.ob<- as.data.frame(out.ob)
out.ob$c <- as.numeric(out.ob$c)
out.ob$z <- as.numeric(out.ob$z)

save(out.pa, file = "C:/Users/Megan Szojka/OneDrive - University of Wyoming/Desktop/Transplants/Scales project/saved_Rdata/mat_pa.RData")
save(out.pb, file = "C:/Users/Megan Szojka/OneDrive - University of Wyoming/Desktop/Transplants/Scales project/saved_Rdata/mat_pb.RData")
save(out.ob, file = "C:/Users/Megan Szojka/OneDrive - University of Wyoming/Desktop/Transplants/Scales project/saved_Rdata/mat_ob.RData")
