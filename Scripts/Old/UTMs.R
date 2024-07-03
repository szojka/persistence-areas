## This script sets up data for persistence area relationship

## Steps
# 1 drop a focal point in one grid on the landscape (take your UTMs and find min   max of x and y, then have r randomly give you 100 places)
# 2 each datasheet will have distance column, where all plots are listed based on   how far away they are from the focal point. Another column will represent        thresholds - plots fall in if their distance is within 25, 25-100 etc.
# 3 do this many times making sure plots have a unique identifier
# 4 last step - link  dataset (run) to seed production

source("Scales project/Scripts/Source_fitness data.R")

#-------------------- Species-Area curve exploration

# see 'Paper1_stats.R' Obj. 3 for these datasets in use

# Species accumulation datasets ####
# set up data like BCI: species are columns, rownames are plots/sites, and values are abundance (for my persistence plots values are seed counts)

# METHOD 1: spaccamum ####
# abiotic persistence sar
pa1 <- dt_p %>%
  dplyr::filter(treatment%in%"A") %>%
  dplyr::select(species,names,seed)%>%
  distinct() 
pa <- pivot_wider(pa1,names_from=species,values_from=seed)# ignore warning
pa <- dplyr::select(pa,-names)
pa <- as.data.frame(pa)
pa$plaere <- as.numeric(as.character(pa$plaere))# na introduced = no prob
pa$brohor <- as.numeric(as.character(pa$brohor))
pa$vulmic <- as.numeric(as.character(pa$vulmic))
pa$miccal <- as.numeric(as.character(pa$miccal))

pa[is.na(pa)] <- 0 # same form as data(BCI)

# biotic persistence sar
pb1 <- dt_p %>%
  dplyr::filter(treatment%in%"B") %>%
  dplyr::select(species,names,seed)%>%
  distinct() 
pb <- pivot_wider(pb1,names_from=species,values_from=seed)
pb <- dplyr::select(pb,-names)
pb <- as.data.frame(pb)
pb$plaere <- as.numeric(as.character(pb$plaere)) 
pb$brohor <- as.numeric(as.character(pb$brohor))
pb$vulmic <- as.numeric(as.character(pb$vulmic))
pb$miccal <- as.numeric(as.character(pb$miccal))

pb[is.na(pb)] <- 0 # same formas data(BCI)

# abiotic occurance sar
oa1 <- dt_o %>%
  dplyr::filter(treatment%in%"A") %>%
  dplyr::select(species,names,ab_cat)%>%
  distinct() 
oa <- pivot_wider(oa1,names_from=species,values_from=ab_cat)
oa <- dplyr::select(oa,-names)
oa <- as.data.frame(oa)
oa$plaere <- as.numeric(as.character(oa$plaere))
oa$brohor <- as.numeric(as.character(oa$brohor))
oa$vulmic <- as.numeric(as.character(oa$vulmic))
oa$miccal <- as.numeric(as.character(oa$miccal))

oa[is.na(oa)] <- 0 # same formas data(BCI)

# biotic occurance sar
ob1 <- dt_o %>%
  dplyr::filter(treatment%in%"B") %>%
  dplyr::select(species,names,ab_cat)%>%
  distinct() 
ob <- pivot_wider(ob1,names_from=species,values_from=ab_cat)
ob <- dplyr::select(ob,-names)
ob <- as.data.frame(ob)
ob$plaere <- as.numeric(as.character(ob$plaere))
ob$brohor <- as.numeric(as.character(ob$brohor))
ob$vulmic <- as.numeric(as.character(ob$vulmic))
ob$miccal <- as.numeric(as.character(ob$miccal))

ob[is.na(ob)] <- 0 # same formas data(BCI)
rm(pa1,pb1,oa1,ob1)
# see 'Paper1_stats.R' for curves

# Species-area relationship datasets ####
# METHOD 2: SSarrhenius ####
# add a unique identifying to each area, this will be equivilant to 'islands' in example data
# dt_o & dt_p

# JUST ALTER MAT_PB MAT_PA MAT_OA AND MAT_OB!! (APPLIES TO OBJ 3 USING METHOD 2 & 3)

# dt_o_mod <- dt_o %>%
#   dplyr::filter(occurance == 'yes')
# dt_map <- dplyr::select(dt_o_mod, area)
# dt_map$area <- round(dt_map$area, digits = 0)
# dt_map <- arrange(dt_map, area)
# class(dt_map$area)
# dt_map$area<- as.factor(as.numeric(dt_map$area))
# # create new column to represent unique areas (isl)
# length(unique(dt_map$area)) #9697
# dt_map <- distinct(dt_map)
# #dt_map$isl <- c(1:9697)
# dt_map <- as.data.frame(dt_map)
# rownames(dt_map) <- dt_map$area
# dt_map$area<- as.numeric(as.character(dt_map$area))
# 
# # now creat equivilant to sipoo
# # occurance biotic
# sarOB2 <- dt_o %>%
#   dplyr::filter(treatment%in%"B") %>%
#   dplyr::select(species,area,ab_cat)%>%
#   distinct() 
# a_ob <- sarOB2$area # to use to filter dt_map later so that areas are identical lengths
# sarOB2 <- arrange(sarOB2, area)
# sarOB2$area <- as.factor(as.numeric(sarOB2$area))
# sarOB2$ab_cat[sarOB2$ab_cat>0] <- 1
# 
# sarOB1 <- as.data.frame(sarOB2)
# sarOB1 <- distinct(sarOB1)
# sarOB1 <- pivot_wider(sarOB1,names_from="species",values_from="ab_cat")
# sarOB1[is.na(sarOB1)] <- 0
# sarOB1 <- as.data.frame(sarOB1)
# sarOB1$area<- as.numeric(as.character(sarOB1$area))
# 
# sarOB <- sarOB1
# rownames(sarOB) <- sarOB$area
# sarOB <- dplyr::select(sarOB, -area)
# 
# sarOB$plaere <- as.numeric(as.character(sarOB$plaere))
# sarOB$brohor <- as.numeric(as.character(sarOB$brohor))
# sarOB$vulmic <- as.numeric(as.character(sarOB$vulmic))
# sarOB$miccal <- as.numeric(as.character(sarOB$miccal))
# # order sarOB from smallest rowSums to largest
# 
# class(sarOB)
# length(sarOB)
# sarOB_mat <- as.matrix(sarOB)
# class(sarOB_mat)
# length(sarOB_mat) #9694
# length(dt_map$area) #9697
# # there is a 3 area difference, find it!
# length(unique(dt_o$area)) #[1] 9697
# # because sarOB is biotic and occurance subset not all scales will be represented. 
# # solution:
# dt_map_ob <- dplyr::filter(dt_map, area%in%a_ob)

# REPEAT FOR A and P combos


# METHOD 3: mmSAR #### 
# mmOB1 <- sarOB1
# mmOB1$s <- rowSums(mmOB1)
# mmOB1 <- dplyr::select(mmOB1, -c(2:5))
# names(mmOB1)[1] <- "a" 
# name = 'Arrhenius 1921'
# data = mmOB1
# mmOB <- list(name, data)
# names(mmOB) <- c("name","data")

#-----------------------------------------------------------------------.
# Response Variable Calculation ####
# the below loop calculates the average # of peristing species as area increases
nrows <- 15*2*nrun
mat <- matrix(data = NA, nrow = nrows, ncol = 4) 
colnames(mat) <- c("run", "scale", "treatment", "sp")

k <- 1 # counter for vector for matrix 'mat' rows

for (r in unique(df$run)) {
  temp <- subset(df, run == r)
  # because scale levels often change between runs
  
  n_temp <- subset(n_sum, run == r)
  # so that averages reflect the correct sum within each run
  
  for (i in unique(df$scale)) {
    temp2 <- subset(temp, scale == i)
    n_temp2 <- subset(n_temp, scale == i)
    
    for (t in levels(df$treatment)) {
      # calculates for A first then B (see arrange in df) as keeping an order is important for later
      
      sp <- # scale must be numeric for this to work
        sum(temp2$persist[temp2$scale <= i & temp2$treatment == t]) /
        sum(n_temp2$n[n_temp2$scale <= i & n_temp2$treatment == t])
      # dived by the sum of n for the given scale & treatment to get an appropriate mean at each group
      
      # fill the matrix with current values
      mat[k,1] <- r
      mat[k,2] <- i
      mat[k,3] <- t
      mat[k,4] <- sp
      
      k <- k + 1
    }
  }
} 

# fills to exactly the data that exists
mat.df <- as.data.frame(mat)
mat.df <- dplyr::filter(mat.df, !sp %in% "NaN")
str(mat.df)

mat.df$treatment <- as.factor(as.character(mat.df$treatment))
mat.df$scale <- as.numeric(as.character(mat.df$scale))
mat.df$sp <- as.numeric(as.character(mat.df$sp))

# combine mat.df & dat
dat <- left_join(dat, mat.df, by = c("run","scale","treatment"))

# rescale sp where 0 = 0% and 4=100%
dat$sp <- rescale(dat$sp, from = c(0,4), to = c(0,100)) 

rm(mat.df, n_sum, temp, n_temp, temp2, n_temp2, dat1, dist.dat, dt1, dt2)

# df.scales ####
# find distances between each of my plot points in "utm2"
y <- pointDistance(utm2, lonlat= FALSE, allpairs = TRUE)
dist.dat <- as.data.frame(t(y))
names <- c(utm$name)# create vector of names from utm

# collapse dataframe (dist.dat) into two columns where column names are third column indicating value 
dat1 <- pivot_longer(dist.dat, cols = c(1:450), names_to = "loc")

# prepare dat1 to join with new_key below
dat1$names <- rep(names, each=450)
names(dat1)[2] <- c("area")
# scale categories for test
dat1$scale <- NA
dat1$scale[dat1$area <= 50] <- 'local'
dat1$scale[dat1$area <= 200 & dat1$area > 50] <- 'regional'
dat1$scale[dat1$area > 200] <- 'global'
table(dat1$scale)
dat1 <- dplyr::select(dat1, -area)

# response variable 'n' is added here
# plot n
pn <- plotlev %>% 
  group_by(treatment, species, contingency) %>%
  add_count() %>%
  distinct()
pn$scale <- "local" 

pn <- left_join(pn, name_key, by = c('tag','treatment','species','seed','grid',"site"))
pn <- dplyr::select(pn, treatment, species, contingency, names, n, scale)

# grid n
gn <- gridlev %>% 
  group_by(treatment, species, contingency) %>%
  add_count() %>%
  distinct()
gn$scale <- "regional" 

gn <- left_join(gn, name_key, by = c('treatment','species','grid',"site"))
gn <- dplyr::select(gn, treatment, species, contingency, names, n, scale)

# site n
sn <- sitelev %>% 
  group_by(treatment, species, contingency) %>%
  add_count() %>%
  distinct()
sn$scale <- "global" 

sn <- left_join(sn, name_key, by = c('treatment','species',"site"))
sn <- dplyr::select(sn, treatment, species, contingency, names, n, scale)

new_key <- rbind(pn,gn,sn)

# join
df.scales <- left_join(dat1, new_key, c("names","scale"))
df.scales <- distinct(df.scales)

# checking all is well with this dataset
str(df.scales)
df.scales$loc <- as.factor(as.character(df.scales$loc))
levels(df.scales$loc) # good: 450 exist
df.scales$treatment <- as.factor(as.character(df.scales$treatment))
levels(df.scales$treatment) # good
df.scales$scale  <- as.factor(as.character(df.scales$scale))
df.scales <- na.omit(df.scales)

rm(y, dat1, dist.dat)

df.scales <- df.scales %>%
  dplyr::select(-names) %>%
  distinct() %>%
  filter(loc%in%"V387") %>%
  dplyr::select(-loc)
