
#-------------------------------------------------------------------------------
# DESCRIPTION: create SAR PAR categories using averaged suitable habitat: 
# diversity, realized-PAR and potential-PAR both with and without neighbors.
# averaged persistence calculated as seeds in v seeds out at plot, grid, site scales.
#-------------------------------------------------------------------------------

# using plotlev, gridlev, sitelev

##############################
# Load packages and functions
##############################

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

#--------------------------------------------
# Plot level
#--------------------------------------------

########################
# Data wrangling
########################

# For poa and pob 
# need to have separate dataframes for just O info and just P info. NOT connected as type.
dt_p <- dt1 %>% # dt1 is plotlev (length 3060) at 
  dplyr::select(-occurrence, -ab_cat) #%>% dplyr::filter(!persistence%in% 'no')
dt_o <- dt1 %>%
  dplyr::select(-persistence, -seed) #%>% dplyr::filter(!occurrence%in%'no')
dt_op <- left_join(dt_o, dt_p, by = c('species','treatment','names','tag', 'grid', 'site'))
dt_op <- dt_op %>% #dplyr::filter(occurrence == 'yes') %>% # removed this because we want the zeros that come from occ = no
  dplyr::select(-ab_cat, -tag, -seed)
dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
dt_op$occurrence <- as.numeric(dt_op$occurrence)
dt_op$persistence[dt_op$persistence=='yes'] <- 1
dt_op$persistence[dt_op$persistence == 'no'] <- 0
dt_op$persistence<- as.numeric(dt_op$persistence)

# if statements: if both occ and persistence = 0, then 0, if occ = 0, then 0, if persist = 0 , then 0, if pers AND occ have 1, then 1
dt_op$realized <- NA
dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                        ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))

dt_op <- dt_op %>% 
  select(-occurrence, -persistence) %>%
  distinct() %>%
  group_by(names,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
  dplyr::mutate(replicates = n()) %>%
  tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species

# NA to zeros
dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
dt_op$species_no <- rowSums(dt_op[c(6:9)])

dt_op$grid <- as.factor(dt_op$grid)
dt_op$site <- as.factor(dt_op$site)
dt_op$treatment <- as.factor(dt_op$treatment)

# good data for POA or POB
  
# for OB PA and PB  
temp <- plotlev # easier coding than yes/no
  temp$persistence <- case_match( temp$persistence, 
                                    'yes' ~ "1",
                                    'no' ~ "0",
                                    .default =  temp$persistence)
  temp$occurrence <- case_match( temp$occurrence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$occurrence)
  temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(names,  species, treatment, persistence, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(names, treatment) %>%
    dplyr::mutate(replicates = n()) %>%
    distinct()
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
  temp_p[is.na(temp_p)] <- 0
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(6:9)])
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(names, species, treatment, occurrence, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(names, treatment) %>%
    dplyr::mutate(replicates = n()) %>%
    distinct()
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(6:9)])
 
  temp_o$type <- as.factor(temp_o$type)
  temp_p$type <- as.factor(temp_p$type)

  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment) %>%
    dplyr::distinct()
  
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()

  
  # this matters for how I want to fit the models. Treatment as predictor or separate??
  # I think type as predictor, with treatment implicit.
  plot_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  plot_par <- plot_par %>%
    dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
    mutate(species_prop = species_no/replicates)
  
  # make type a factor
  plot_par$type <- as.factor(plot_par$type)
  plot_par$scale <- c("plot")

#--------------------------------------------
# Grid level
#--------------------------------------------

########################
# Data wrangling
########################

  # Create 'realized' data for poa and pob 
  
  # need to have separate columns for just O info and just P info. NOT connected as type.
  dt_op <-  gridlev %>% dplyr::select(-grid_n, -grid_seed, -contingency)
  
  dt_op$occurrence <- as.character(dt_op$occurrence)
  dt_op$persistence <- as.character(dt_op$persistence)
  
  dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
  dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
  dt_op$occurrence <- as.numeric(dt_op$occurrence)
  dt_op$persistence[dt_op$persistence=='yes'] <- 1
  dt_op$persistence[dt_op$persistence == 'no'] <- 0
  dt_op$persistence<- as.numeric(dt_op$persistence)
  
  # if statements: if both occ and persistence = 0, then 0, if occ = 0, then 0, if persist = 0 , then 0, if pers AND occ have 1, then 1
  dt_op$realized <- NA
  dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                           ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))
  
  dt_op <- dt_op %>% 
    select(-occurrence, -persistence) %>%
    distinct() %>%
    group_by(grid,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
    dplyr::mutate(replicates = n()) %>%
    tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species
  
  # NA to zeros
  dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
  dt_op$species_no <- rowSums(dt_op[c(5:8)])
  
  dt_op$grid <- as.factor(dt_op$grid)
  dt_op$site <- as.factor(dt_op$site)
  dt_op$treatment <- as.factor(dt_op$treatment)
  # good data for POA or POB
  
 # DATA "POTENTIAL" FOR PA PB
  
  temp <- gridlev # easier coding than yes/no
  temp$persistence <- case_match( temp$persistence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$persistence)
  temp$occurrence <- case_match( temp$occurrence, 
                                 'yes' ~ "1",
                                 'no' ~ "0",
                                 .default =  temp$occurrence)
  temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(species, treatment, persistence, grid, site) %>%
    dplyr::distinct()  %>%
    dplyr::group_by(grid, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
  temp_p[is.na(temp_p)] <- 0 # this is fine, replicates have already been calculated, turn to 0 just to allow rowSum function to work.
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(5:8)]) # only places where there were NAs in resevoir grid 26.
  
  # Check that things are correct is that potential should always be equal or higher than realized
  dt_op %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no)) 
  # 1 A            35
  # 2 B            12
  temp_p %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no))
  # 1 A            37
  # 2 B            12
  # check good.
  
  # DATA "DIVERSITY" FOR OB
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(species, treatment, occurrence, grid, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(grid, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0 
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(5:8)])

  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  # this matters for how I want to fit the models. Treatment as predictor or separate??
  # I think type as predictor, with treatment implicit.
  grid_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  grid_par <- grid_par %>%
    dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
    mutate(species_prop = species_no/replicates) 
  
  # make type a factor
  grid_par$type <- as.factor(grid_par$type)
  grid_par$scale <- c("grid")
  
#--------------------------------------------
# Site level
#--------------------------------------------

########################
# Data wrangling
########################

  # For poa and pob 
  # need to have separate columns for just O info and just P info. NOT connected as type.
  dt_op <- sitelev %>% dplyr::select(-site_n, -site_seed, -contingency) 

  dt_op$occurrence <- as.character(dt_op$occurrence)
  dt_op$persistence <- as.character(dt_op$persistence)
  
  dt_op$occurrence[dt_op$occurrence=='yes'] <- 1
  dt_op$occurrence[dt_op$occurrence == 'no'] <- 0
  dt_op$occurrence <- as.numeric(dt_op$occurrence)
  dt_op$persistence[dt_op$persistence=='yes'] <- 1
  dt_op$persistence[dt_op$persistence == 'no'] <- 0
  dt_op$persistence<- as.numeric(dt_op$persistence)
  
  # if statements: if both occ and persistence = 0, then 0, if occ = 0, then 0, if persist = 0 , then 0, if pers AND occ have 1, then 1
  dt_op$realized <- NA
  dt_op$realized <- ifelse(dt_op$occurrence == 0 | dt_op$persistence == 0, 0, 
                           ifelse(dt_op$occurrence == 1 & dt_op$persistence == 1, 1, NA))
  
  dt_op <- dt_op %>% 
    dplyr::select(-occurrence, -persistence) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site,treatment) %>%# so now my replicates are count species at each name and treatment that are both persisting and occurring
    dplyr::mutate(replicates = n()) %>%
    tidyr::pivot_wider(names_from = species, values_from = realized) # pivot to count persisting & occurring species
   
  # NA to zeros
  dt_op[is.na(dt_op)] <- 0 # must turn these to 0 for row sums, but already calculated replicates so no problem
  dt_op$species_no <- rowSums(dt_op[c(4:7)])
  
  dt_op$site <- as.factor(dt_op$site)
  dt_op$treatment <- as.factor(dt_op$treatment)
  # good data for POA or POB
  
  # for PA PB OB
  temp <- sitelev %>% # easier coding than yes/no
    dplyr::select(-contingency, -site_n, -site_seed)
  temp$persistence <- case_match( temp$persistence, 
                                  'yes' ~ "1",
                                  'no' ~ "0",
                                  .default =  temp$persistence)
  temp$occurrence <- case_match( temp$occurrence, 
                                 'yes' ~ "1",
                                 #'no' ~ "0",
                                 .default =  temp$occurrence)
  temp$persistence <- as.numeric(temp$persistence)
  temp$occurrence <- as.numeric(temp$occurrence)
  
  temp_p <- temp %>% # temporary data frame to manipulate persistence (yes or no)
    dplyr::select(species, treatment, persistence, site) %>%
    dplyr::distinct()  %>%
    dplyr::group_by(site, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_p <- pivot_wider(temp_p, names_from = "species", values_from = "persistence")
  temp_p[is.na(temp_p)] <- 0
  temp_p$type <- "P"
  temp_p$species_no <- rowSums(temp_p[c(4:7)])
  
  temp_o <- temp %>% # temporary data frame to manipulate occurrence (yes or no)
    dplyr::select(species, treatment, occurrence, site) %>%
    dplyr::distinct() %>%
    dplyr::group_by(site, treatment) %>%
    dplyr::mutate(replicates = n())
  temp_o <- pivot_wider(temp_o, names_from = "species", values_from = "occurrence")
  temp_o[is.na(temp_o)] <- 0
  temp_o$type <- "O"
  temp_o$species_no <- rowSums(temp_o[c(4:7)])
  
  # combine all types together, using persistence and occurrence but getting more specific.
  pob <- dt_op %>% 
    dplyr::filter(treatment %in% c('B')) %>%
    dplyr::mutate(type = "POB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  poa <- dt_op %>% 
    dplyr::filter(treatment %in% c('A')) %>%
    dplyr::mutate(type = "POA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  ob <- temp_o %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('O')) %>%
    dplyr::mutate(type = "OB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pb <- temp_p %>% 
    dplyr::filter(treatment %in% c('B') & type %in% c('P')) %>%
    dplyr::mutate(type = "PB") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  pa <- temp_p %>% 
    dplyr::filter(treatment %in% c('A') & type %in% c('P')) %>%
    dplyr::mutate(type = "PA") %>%
    dplyr::ungroup() %>%
    dplyr::select(-treatment)%>%
    dplyr::distinct()
  
  # this matters for how I want to fit the models. Treatment as predictor or separate??
  # I think type as predictor, with treatment implicit.
  site_par <- rbind(pob,ob,pb,poa,pa) # all these are types.
  
  # make response variable proportions and remove species
  site_par <- site_par %>%
  dplyr::select(-plaere,-miccal,-brohor,-vulmic) %>%
  dplyr::mutate(species_prop = species_no/replicates) 
  
  # make type a factor
  site_par$type <- as.factor(site_par$type)
  site_par$scale <- c("site")

  # Check that things are correct is that potential should always be equal or higher than realized
  dt_op %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no)) 
  # 1 A            14
  # 2 B             4
  temp_p %>%
    dplyr::group_by(treatment) %>%
    dplyr::summarize(tot = sum(species_no))
  # 1 A            14
  # 2 B             4
  # check good, exactly the same.
  
################################################
# BRING ALL RAW DATAFRAMES TOGETHER FOR PLOTTING

grid_par$names <- NA
site_par$names <- NA
site_par$grid <- NA

names(plot_par)
names(grid_par)
names(site_par)
par_all_discrete <- rbind(plot_par, grid_par, site_par)

# Change variable types / names for all dataframes
par_all_discrete$scale <- as.factor(par_all_discrete$scale)

# order scales plot grid site
par_all_discrete$scale <- factor(par_all_discrete$scale, 
                                 levels = c('plot', 'grid', 'site'))

# rename ob pa poa pob
par_all_discrete$type <- case_match(par_all_discrete$type, 
                                'OB' ~ "Diversity (SAR)",
                                'PA' ~ "Potential \n without neighbors",
                                'POA' ~ "Realized \n without neighbors",
                                'PB' ~ "Potential \n with neighbors",
                                'POB' ~ "Realized \n with neighbors",
                                .default =  par_all_discrete$type)

par_all_discrete$type<- as.factor(par_all_discrete$type)
par_all_discrete$names<- as.factor(par_all_discrete$names)

############################
# summary stats and checks
############################
par_all_discrete <- par_all_discrete %>%
  dplyr::group_by(type,scale) %>%
  dplyr::mutate(mean = mean(species_prop)) 
ggplot(par_all_discrete) + 
 geom_histogram(aes(x=species_prop)) +
  facet_grid(rows = vars(scale), cols = vars(type)) + # realized are both zero inflated or more poisson like
  geom_vline(aes(xintercept = mean), color = "red")

par_all_discrete %>%
  dplyr::group_by(type, scale) %>%
  dplyr::mutate(min = min(species_prop)) %>%
  dplyr::mutate(max = max(species_prop)) %>%
  dplyr::mutate(std = sd(species_prop)) %>%
  dplyr::mutate(N = n()) %>%
  dplyr::select(type, scale, min, max, std, mean, N) %>%
  distinct()

# type                             scale   min   max   std   mean
# <fct>                            <fct> <dbl> <dbl> <dbl>  <dbl>
# 1 "Realized \n with neighbors"     plot    0    1    0.171 0.111   443
# 2 "Diversity (SAR)"                plot    0    1    0.305 0.561   443
# 3 "Potential \n with neighbors"    plot    0    1    0.217 0.190   443
# 4 "Realized \n without neighbors"  plot    0    1    0.247 0.206   422
# 5 "Potential \n without neighbors" plot    0    1    0.269 0.334   422
# 6 "Realized \n with neighbors"     grid    0    0.5  0.224 0.181    18
# 7 "Diversity (SAR)"                grid    0.5  1    0.137 0.944    18
# 8 "Potential \n with neighbors"    grid    0    0.5  0.224 0.181    18
# 9 "Realized \n without neighbors"  grid    0    0.75 0.243 0.5      18
# 10 "Potential \n without neighbors" grid    0    1    0.256 0.528    18
# 11 "Realized \n with neighbors"     site    0    0.5  0.258 0.167     6
# 12 "Diversity (SAR)"                site    1    1    0     1         6
# 13 "Potential \n with neighbors"    site    0    0.5  0.258 0.167     6
# 14 "Realized \n without neighbors"  site    0    0.75 0.303 0.583     6
# 15 "Potential \n without neighbors" site    0    0.75 0.303 0.583     6
  
#----------------------------------------------------------------------------
# FIRST ATTEMPT BELOW WAS MODELLING TYPE AS FACTOR RESPONSE.

###############
# Fit models 
###############

############################
# ABIOTIC & BIOTIC TOGETHER

#------------
# Plot
# 
m.p <- glmmTMB::glmmTMB(species_prop ~ type + (1|names) + (1|site:grid) + (1|site),
                        data = dplyr::filter(par_all_discrete, scale %in% 'plot'),
                        family = binomial(),
                        weights = replicates)

# Check model fit
testZeroInflation(m.p)
testDispersion(m.p)
plot(fitted(m.p), residuals(m.p))
hist(residuals(m.p))
loc_sim_ouput <- simulateResiduals(m.p)
plot(loc_sim_ouput)
testOutliers( # DHARMa:testOutliers with type = binomial may have inflated Type I error rates for integer-valued distributions. To get a more exact result, it is recommended to re-run testOutliers with type = 'bootstrap'. See ?testOutliers for details
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.p)

Anova(m.p, type = 3)
#                 Chisq Df Pr(>Chisq)
# (Intercept)    3.5764  1    0.05861 .  
# type        888.0722  4    < 2e-16 ***

v <- ggpredict(m.p, terms = c("type"), 
               type = "fe", allow.new.levels=TRUE); plot(v)

# add pairwise comparison table
em <- emmeans(m.p, ~type, type = "response") # specify green_index_scaled values
em

p.contrast <- as.data.frame(pairs(em))
pairs(em)

#------------
# Grid

m.g <- glmmTMB::glmmTMB(species_prop ~ type + (1|site:grid) + (1|site),
                        data = dplyr::filter(par_all_discrete, scale %in% 'grid'),
                        family = binomial(),
                        weights = replicates)


# Check model fit
testZeroInflation(m.g)
testDispersion(m.g)
plot(fitted(m.g), residuals(m.g))
hist(residuals(m.g))
loc_sim_ouput <- simulateResiduals(m.g)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) # looks great,


summary(m.g)

Anova(m.g, type = 3)
#               Chisq Df Pr(>Chisq)
#   (Intercept) 53.836  1   2.18e-13 ***
#   type        82.552  2  < 2.2e-16 ***

v <- ggpredict(m.g, terms = c("type"), 
               type = "fe", allow.new.levels=TRUE); plot(v)

# add pairwise comparison table
em <- emmeans(m.g, ~type, type = "response") # specify green_index_scaled values
em

g.contrast <- as.data.frame(pairs(em))
pairs(em)

#------------
# Site

m.s <- glmmTMB::glmmTMB(species_prop ~ type + (1|site),
                        data = dplyr::filter(par_all_discrete, scale %in% 'site'),
                        family = binomial(),
                        weights = replicates)

# Check model fit
testZeroInflation(m.s)
testDispersion(m.s)
plot(fitted(m.s), residuals(m.s))
hist(residuals(m.s))
loc_sim_ouput <- simulateResiduals(m.s)
plot(loc_sim_ouput)
testOutliers(
  loc_sim_ouput,
  alternative = c("two.sided"),
  margin = c("both"),
  type = c("bootstrap"),
  nBoot = 100,
  plot = T
) 

summary(m.s)

Anova(m.s, type = 3)
#               Chisq Df Pr(>Chisq)
# (Intercept)  0.000  1   0.998978   
# type        16.064  4   0.002934 **

v <- ggpredict(m.s, terms = c("type"), 
               type = "fe", allow.new.levels=TRUE); plot(v) 
# all those at 100 must have variances removed. Realized does not!

# add pairwise comparison table
em <- emmeans(m.s, ~type, type = "response") # specify green_index_scaled values
em

s.contrast <- as.data.frame(pairs(em))
pairs(em)

##############################################
# PULL PREDICTIONS TOGETHER WITH GGPREDICT
##############################################

vis2p <- ggpredict(m.p, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis2p$group <- 'plot'

vis2g <- ggpredict(m.g,
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis2g$group <- 'grid'

vis2s <- ggpredict(m.s, 
                   terms = c("type"), 
                   type = "fe", allow.new.levels=TRUE)
vis2s$group <- 'site'
# take confs to 1, as no error can really be estimated when all values are 1.
vis2s$conf.low[vis2s$x %in% 'Diversity (SAR)'] <- 1 

vis2 <- rbind(vis2p,vis2g,vis2s)
vis2$group <- factor(vis2$group, c("plot","grid","site")) # order group into increasing scales

names(vis2) # all correct classes.

vis2 <- as.data.frame(vis2)
# add scenario, kind, treatment, and scale_type

vis2$treatment <- NA
vis2$treatment[vis2$x %in% 
                 c("Potential \n with neighbors","Diversity (SAR)","Realized \n with neighbors")] <- "With neighbors"
vis2$treatment[vis2$x %in% 
                 c("Potential \n without neighbors","Realized \n without neighbors")] <- "Without neighbors"

vis2$kind <- NA
vis2$kind[vis2$x %in% 
            c("Potential \n with neighbors","Potential \n without neighbors")] <- "Potential (PAR)"

vis2$kind[vis2$x %in% 
            c("Realized \n with neighbors","Realized \n without neighbors")] <- "Realized (PAR)"

vis2$kind[vis2$x %in% 
            c("Diversity (SAR)")] <- "Diversity (SAR)"

vis2$scenario <- NA
vis2$scenario <- "Data"

vis2$scale_type <- NA
vis2$scale_type <- "Averaged"

names(vis2)[1] <- "type"
names(vis2)[6] <- "scale"


#############################
# TABLE OUTPUT
#############################

# using gt
#library(gt) # loaded in source

# https://gt.albert-rapp.de/getting_started.html
# necessary for saving but shouldn't have to use more than once.
#Sys.setenv(CHROMOTE_CHROME = "C:/Users/Megan Szojka/AppData/Local/Google/Chrome/Application/chrome.exe")
#chromote::find_chrome()

# vis2, organize by scale, type, predicted, confidence intervals "[{conf.low}, {conf.high}]]
# discard 'treatment' and 'scenario' and 'kind'

paste(round(vis2$conf.low,2), round(vis2$conf.high,1), sep = ", ")

tab_avg_PAR_estimates <- vis2 |>
  dplyr::mutate(conf.int =paste(round(conf.low,2), round(conf.high,1), sep = ", ")) |>
  dplyr::mutate(predicted = round(predicted,2)) |>
 # dplyr::select(-kind, -scenario, -treatment, -std.error, -scale_type, -conf.low, -conf.high) |>
  dplyr::select(type, scale, predicted, conf.int) |>
  gt() |>
  tab_header( title = "",
    subtitle = "Table S12. Model estimates of persistence-area relationships (realized-PAR and potential-PAR) and species-area relationships (Diversity-SAR) using averaged scaling of suitable habitat. Separate models were fit for each scale.")  |>
  opt_align_table_header(align = "left") |>
  cols_label(
    type = 'Data type',
    predicted = 'Estimated proportion',
    conf.int = "95% CI",
    scale = 'Scale level') |>
  cols_align(
    align = 'right', 
    columns = where(is.numeric)) |> 
  cols_align(
    align = 'left', 
    columns = where(is.factor))
  #as_latex() # exports code but don't know how to work it right in latex
  
 tab_avg_PAR_estimates |>
  gtsave(paste0(here::here(),"/Tables/12tab_avg_PAR_estimates.pdf")) 

 ######################
 # 2. ANOVA 
 
 # three models one for each scale
 #m.p, m.g, m.s
 
 a1 <- Anova(m.p, type = 3) 
 a2 <- Anova(m.g, type = 3)
 a3 <- Anova(m.s, type = 3)
 
 anova_avg <- data.frame(Chi.squared = round(c(a1$Chisq[1],a1$Chisq[2],a2$Chisq[1],a2$Chisq[2],a3$Chisq[1],a3$Chisq[2]),3),
                         scale = c('plot','plot','grid','grid','site','site'),
                         Df =c(a1$Df[1],a1$Df[2],a2$Df[1],a2$Df[2],a3$Df[1],a3$Df[2]),
                         P_value =  round(c(a1$`Pr(>Chisq)`[1],a1$`Pr(>Chisq)`[2],a2$`Pr(>Chisq)`[1],a2$`Pr(>Chisq)`[2],a3$`Pr(>Chisq)`[1],a3$`Pr(>Chisq)`[2]),3),
                         predictor = c(rep(c("Intercept","Data type"), times = 3)))
 
 tab_avg_anova <- anova_avg |>
   dplyr::select(scale, predictor, Chi.squared, Df, P_value) |>
   gt() |>
   tab_header( title = "",
               subtitle = "Table S13. ANOVA outputs of predictor 'data type', (which includes persistence-area relationships (realized-PAR and potential-PAR) both with and without neighbors, and species-area relationships (Diversity-SAR)) using averaged scaling of suitable habitat. Separate models were fit for each scale.")  |>
   opt_align_table_header(align = "left") |>
   cols_label(
     predictor = 'Predictor',
     scale = 'Scale level',
     P_value = "P-value") |>
   cols_align(
     align = 'right', 
     columns = where(is.numeric)) |> 
   cols_align(
     align = 'left', 
     columns = where(is.factor))
 
 tab_avg_anova |>
  gtsave(paste0(here::here(),"/Tables/13tab_avg_anova.pdf")) 
 
 ########################
 # 3. pairwise contrasts
 
 contrastz <- rbind(p.contrast, g.contrast, s.contrast)
 contrastz$scale <- NA
 contrastz$scale <- c(rep('plot',times = 10),rep('grid',times = 10), rep('site',times = 10))
 
 levels(contrastz$contrast)
 # remove \n from the level contrasts
 contrastz$contrast <- case_match(contrastz$contrast,  
                                  "Diversity (SAR) / Potential \n with neighbors"  ~ "SAR / Potential-PAR with neighbors",            
                                  "Diversity (SAR) / Potential \n without neighbors"  ~ "SAR / Potential-PAR without neighbors",       
                                  "Diversity (SAR) / Realized \n with neighbors"  ~  "SAR / Realized-PAR with neighbors",             
                                  "Diversity (SAR) / Realized \n without neighbors" ~ "SAR / Realized-PAR without neighbors",
                                  "Potential \n with neighbors / Potential \n without neighbors"  ~ "Potential-PAR with neighbors / Potential-PAR without neighbors",
                                  "Potential \n with neighbors / Realized \n with neighbors"  ~ "Potential-PAR with neighbors / Realized-PAR with neighbors",     
                                  "Potential \n with neighbors / Realized \n without neighbors" ~ "Potential-PAR with neighbors / Realized-PAR without neighbors" ,  
                                  "Potential \n without neighbors / Realized \n with neighbors"  ~  "Potential-PAR without neighbors / Realized-PAR with neighbors", 
                                  "Potential \n without neighbors / Realized \n without neighbors" ~ "Potential-PAR without neighbors / Realized-PAR without neighbors",
                                  "Realized \n with neighbors / Realized \n without neighbors"  ~"Realized-PAR with neighbors / Realized-PAR without neighbors" ) 
 
 tab_avgPAR_contrasts <- contrastz |>
   dplyr::select(scale, contrast, z.ratio, p.value) |>
   # dplyr::mutate(odds.ratio = round(odds.ratio,3)) |>
   dplyr::mutate(z.ratio = round(z.ratio,3)) |>
   dplyr::mutate(p.value = round(p.value,3)) |>
   gt() |>
   tab_header( title = "",
               subtitle = "Table S14. Contrasts between each PAR and SAR from averaged method of scaling suitable habitat. Each scale was fit with a separate model.")  |>
   opt_align_table_header(align = "left") |>
   cols_label(
     p.value = "P-value",
     contrast = 'Pairwise comparison',
     scale = "Scale level",
     z.ratio = "z-ratio") |>
   cols_align(
     align = 'right', 
     columns = where(is.numeric)) |> 
   cols_align(
     align = 'left', 
     columns = where(is.factor))
 
tab_avgPAR_contrasts |>
  gtsave(paste0(here::here(),"/Tables/14tab_avg_PAR_contrasts.pdf")) 
 
 # EXTRA CODE THAT IS MORE ELEGANT THAN MULTISTEP ABOVE
 # tt <- sitelev %>% # this works for averaged way, delete
 #   dplyr::filter(treatment %in% 'B') %>%
 #   dplyr::arrange(site) %>%
 #   dplyr::group_by(site) %>%
 #   dplyr::mutate(replicates = n()) %>%
 #   dplyr::select(-contingency, -site_n, -site_seed) %>%
 #   dplyr::mutate(realized_site = ifelse(occurrence %in% 'yes' & persistence %in% 'yes', 1, 0)) %>%
 #   dplyr::mutate(potential_site = ifelse(persistence %in% 'yes', 1, 0)) %>%
 #   dplyr::mutate(occurrence = ifelse(occurrence %in% 'yes', 1, 0)) %>%
 #   dplyr::mutate(realized_site_sum = sum(realized_site)) %>% # add up all the species for each site
 #   dplyr::mutate(potential_site_sum = sum(potential_site)) %>% # add up all the species for each site
 #   dplyr::select(-species) %>%
 #   distinct() %>% # remove all the species level information, now at community level
 #   dplyr::mutate(realized_site_prop = realized_site_sum/replicates) %>%
 #   dplyr::mutate(potential_site_prop = potential_site_sum/replicates) %>%
 #   dplyr::select(site, treatment, occurrence, persistence, realized_site, realized_site_prop, potential_site, potential_site_prop, replicates) %>%
 #   print( n = 24)
 #   view(tt)