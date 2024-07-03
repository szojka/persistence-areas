


# DESCRIPTION: 



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
library(ggridges)

source("Scripts/Source - MAIN fitnessdata.R")


## show variability of seed production across scales

# smallest scale = individuals: graph every seed observation by frequency
# using seeds_clean instead of plotlev as it has more obs due to NAs introduced by abundance data
#f1 <-  
  ggplot(seeds_clean, aes(x=seed)) +
  geom_density(aes(fill = treatment), position = "stack",alpha=0.4,size=0.8) +
  scale_fill_manual( values = c("aquamarine","blue3")) +
  theme_classic() +
  theme(legend.position = 'none') +
  geom_vline(aes(xintercept = mean(seed)),
             linetype = "dashed", color = "red") +
  labs(x = "Per capita seed production", y="") +
  xlim(-4,30) 
# goes to seed = 108 so I cut to 30 for viz

# grids
grids <- as.numeric(levels(factor(seeds_clean$grid)))
g <- data.frame()
y <- list()
x <- list()
# working on nested loop version for learning - runs :)
# goes to grid 1, takes all treatment A, takes mean of each species for treatment A, goes to grid 1, takes all treatment B, takes mean of each species ofr treatment B, goes to grid 2...
for(l in grids) {
  treat <- levels(factor(seeds_clean$treatment[seeds_clean$grid == l]))
  for(t in treat) {
    spp <- levels(factor(seeds_clean$species[seeds_clean$treatment == t]))
    for(s in spp) {
      filter(seeds_clean, grid == l, treatment == t, species == s) -> x
      y <- mean(x$seed)
      g <- rbind(g, y)
      print(g)
    }
  }
}
names(g)[1] <- c("fitness")
g$fitness <- round(g$fitness)
head(g)

t_grid <- rep(LETTERS[1:2], each = 4, times = 18)
length(t_grid)
t_grid <- as.data.frame(t_grid)

g <- merge(t_grid, g, by = 0)

#f2 <- 
  ggplot(g, aes(x=fitness)) +
    geom_density(aes(fill = t_grid), position = "stack",alpha=0.4,size=0.8) +
    scale_fill_manual( values = c("aquamarine","blue3")) +
    theme_classic() +
    theme(legend.position = 'none') +
    geom_vline(aes(xintercept = mean(fitness)),
               linetype = "dashed", color = "red") +
    labs(x = "Per capita seed production", y="") +
    xlim(-4,30) 

# sites
seeds_clean$site <- as.factor(as.character(seeds_clean$site))
sites <- levels(factor(seeds_clean$site))  
k <- data.frame() 
x <- list() 
y = list() # something is wrong in the loop as the number of replicates I have doesn't add up to what I should have
for(l in sites) {
  treat <- levels(factor(seeds_clean$treatment[seeds_clean$site == l]))
  for(t in treat) {
    spp <- levels(factor(seeds_clean$species[seeds_clean$treatment == t]))
    for(s in spp) {
      filter(seeds_clean, site == l, treatment == t, species == s) -> x
      y <- mean(x$seed)
      k <- rbind(k, y)
      print(k) # ouput has 48 unique values as it should ! :)
    }
  }
}
names(k)[1] <- c("fitness")
k$fitness <- round(k$fitness)
head(k)

t_site <- rep(LETTERS[1:2], each = 4, times = 6)
length(t_site)
t_site <- as.data.frame(t_site)

k <- merge(t_site, k, by = 0)

#f3 <- 
ggplot(k, aes(x=fitness)) +
  geom_density(aes(fill = t_site), position = "stack",alpha=0.4,size=0.8) +
  scale_fill_manual(values = c("aquamarine","blue3")) +
  theme_classic() +
  theme(legend.position = 'none') +
  geom_vline(aes(xintercept = mean(fitness)),
             linetype = "dashed", color = "red") +
  labs(x = "Per capita seed production", y="") +
  xlim(-4,30) 

# regional scale
r <- seeds_clean %>%
  dplyr::group_by(treatment, species) %>%
  dplyr::summarize(x = mean(seed)) %>%
  print()
r <- as.data.frame(r)
r$x <- round(r$x)
head(r) # n = 8 correct :)

f4 <- ggplot(r, aes(x, fill = treatment)) +
  geom_density(position = "stack") +
  theme_classic() +
  geom_vline(aes(xintercept = mean(x)),
             linetype = "dashed", color = "red") +
  labs(x = "seeds by region")+
  xlim(-1,30)

# bring together
s <- dplyr::select(seeds_clean, treatment, seed)
names(s)[2] <- "fitness"
s$scale <- "1m2"
g <- dplyr::select(g, t_grid, fitness)
names(g)[1] <- 'treatment'
g$scale <- "25m2"
k <- dplyr::select(k, t_site, fitness)
names(k)[1] <- 'treatment'
k$scale <- '100m2'
r <- dplyr::select(r, treatment,x)
names(r)[2] <- 'fitness'
r$scale <- '10000m2'

seeds_scale <- rbind(g,s)
rm(s,g,k,r)
seeds_scale$fit_log <- log(seeds_scale$fitness + 1)
seeds_scale$scale <- factor(seeds_scale$scale, levels = c("25m2",'1m2'))# order scale levels

png("Figures/fitness_dist.png", height = 5, width = 6, "in", res = 300)
ggplot(seeds_scale, aes(x = fit_log, y=scale))+
  geom_density_ridges(scale = 1.5, alpha = 0.4, fill = 'blue3') + 
  scale_y_discrete(expand = c(0, 0)) +     # will generally have to set the `expand` option
  scale_x_continuous(expand = c(0, 0)) +   # for both axes to remove unneeded padding
  coord_cartesian(clip = "off") + # to avoid clipping of the very top of the top ridgeline
  theme_classic()+
  labs(x='Log fitness',y ='Spatial scale')+
  facet_wrap(~treatment) 
dev.off()
# mean has not changed across scales which is a good check 
# the districution has changed

# summary stats:
m <- filter(seeds_scale, scale%in%'1m2'& treatment %in% 'B')
mean(m$fitness)
length(m$fitness)
sd(m$fitness/sqrt(1650))

m <- filter(seeds_scale, scale%in%'25m2'& treatment %in% 'B')
mean(m$fitness)
length(m$fitness)
sd(m$fitness/sqrt(72))

m <- filter(seeds_scale, scale%in%'1m2'& treatment %in% 'A')
mean(m$fitness)
length(m$fitness)
sd(m$fitness/sqrt(1478))

m <- filter(seeds_scale, scale%in%'25m2'& treatment %in% 'A')
mean(m$fitness)
length(m$fitness)
sd(m$fitness/sqrt(72))
