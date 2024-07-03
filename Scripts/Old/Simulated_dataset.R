## This is the pre-data sheet for code to get ready for project ##

# Variables: site (18),  tag (3600), species(1 of 4), treatment (biotic or abiotic), seed, plot (which square - should be 8 individuals per plot, and could be #1-25 in the grid)

# Catories = species -> LETTERS(1:4), treatment -> letters(1:2), Grid -> 1:18
# 25 observations in each grid

# Numerical = seed -> between 0 and 10, distance -> pairwise distances bw grids (maybe)

## Attempt one - seed #, species, treatment - one grid ##

grid.df <- data.frame(species = rep(letters[23:26], each = 50),
           treatment = rep(LETTERS[1:2], times = 2),
           seed = rpois(n = 200, lambda = 3))
summary(grid.df) # max for seed is 8 min is 0, showed that rpois worked

# Two - add grids - each of above must be replicated within grids, so times = 18

df <- data.frame(species = rep(LETTERS[23:26], each = 900),
           treatment = rep(LETTERS[1:2], times = 2),
           seed = rpois(n = 3600, lambda = 3),
           grid = rep(1:18, times = 1))
summary(df) # 200 obs in each grid, looks good!
# for now assume grid # indicates productivity!

library(tidyverse)
# what we would expect with just statistical variation of seeds:
# regionally:
ggplot(df, aes(x = species,y = seed)) +
  geom_boxplot() +
  theme_classic()
# within grids:
ggplot(df, aes(x = species,y = seed)) +
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~grid)

# add a column that adds mean and variation for each grid / species
df <- df %>%
  group_by(grid, species) %>%
  mutate(pop.fit = mean(seed)) %>%
  mutate(pop.var = var(seed))
  
# plot curves of species based on mean and var over x =?
ggplot(df, aes(x=grid,y=seed, color=species)) +
  geom_point() +
  facet_wrap(~species)

# r trunc norm





