# Autocorrelation ####
#   C) ?qualitative? Is the grain of heterogeneity different > than or < than average species dispersal distances? 

# Moran's I
# Question is how related the values of a variable are based on the locations where they were measured?
# REDO WITH SOME VALUES FROM MATRIX

# Steps:

# 1) use buffer's for value (mean_ndvi in df) for each site using only year 2016
locs <- read_csv("McLaughlin site info/szojka_locs_2020.csv")
names(locs)[2] <- "grid"
locs <- dplyr::select(locs, grid, latitude, longitude)

p <- left_join(ndvi_mean, locs, by = "grid")

# 2) create inverse distance matrix

library(ape)

pt_dists <- as.matrix(dist(cbind(p$longitude, p$latitude)))
pt_dists_inv <- 1/pt_dists
diag(pt_dists_inv) <- 0

# 3) calculate Moran's I using ape package

spatial_output <- Moran.I(p$mean_ndvi, pt_dists_inv)
spatial_output
# 
# $observed
# [1] 0.7366904
# 
# $expected
# [1] -0.05882353
# 
# $sd
# [1] 0.1891152
# 
# $p.value
# [1] 2.593502e-05
#
# reject Null Hyp that zero spatial autocorrelation is present

# with dispersal distances...

# Supplemental ####
# relationship between seed production and species
class(seeds_clean$grid) # numeric
seeds_clean$grid <- as.factor(as.numeric(seeds_clean$grid))

m_sup <- lm(seed ~ species + grid, data = seeds_clean)
summary(m_sup)
Anova(m_sup, type = 3)
# grids, and species, produce very different amount of seeds

  