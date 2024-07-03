# Obj 1 ####
# How frequently does realized niches differ from fundamental niches, and what is the role of biotic interactions in changing these?

# Test 1) thinking only of communities (==grids, n=18) …
# use Chi-test of independence (Anova works on means and variance, I have observed n & expected n)

# does the frequency of community statuses differ between biotic and abiotic treatments?

## METHOD 1 - Muli reg test ####
# using Rachel's fake df.R example
# plot = names
# species = species
# type = contingency
# code = names + species
# random = names + species + contingency (I call it unique)

#### local level ####
# random effect is grid|site at local scale
# random effect is site|1 at grid level
# random effect is nothing at site level

# local abiotic stats
lm<-glmmTMB(data ~ contingency*species + (site|grid), data=dla, family="binomial") # or should it be...
#lm<-glmmTMB(data ~ contingency + (1|species), data=dfa,family="binomial")
Anova(lm)
vis<-visreg(lm,xvar="species",scale="response")
summary(lm)
vis$fit$visregFit
sum(vis$fit$visregFit) # sums to 1 ish :)

# local biotic stats
lm<-glmmTMB(data ~ contingency*species + (site|grid), data=dlb, family="binomial") # or should it be...
#lm<-glmmTMB(data ~ contingency + (1|species), data=dfa,family="binomial")
Anova(lm)
vis<-visreg(lm,xvar="contingency",scale="response")

vis$fit$visregFit
sum(vis$fit$visregFit) # sums to 1 ish :)

# grid abiotic stats
lm2<-glmmTMB(data ~ contingency*species + (1|site), data=dga, family="binomial") 
#non convergence
#taking away random effect lets it run
# trouble shooting with Control function
lm2.1 <- update(lm2, control=glmmTMBControl(optimizer=optim,
                                            optArgs=list(method="BFGS")))
lm2.1$fit$convergence # 0, meaning convergence was successful (Rdoc)
Anova(lm2.1) 
vis<-visreg(lm2,xvar="contingency",scale="response")

vis$fit$visregFit
sum(vis$fit$visregFit)

# site abiotic stats
#non convergence
lm3<-glmmTMB(data ~ contingency*species, data=dsa, family="binomial")
# trouble shooting with Control function
lm3.1 <- update(lm3, control=glmmTMBControl(optimizer=optim,
                                            optArgs=list(method="BFGS"))) # warning for extremely small eigen values detected
lm3.1$fit$convergence # 0, meaning it converged
Anova(lm3.1)
vis<-visreg(lm3.1,xvar="contingency",scale="response") 
##### PROBLEM #
# model says that only SS_y is predicted, though some are ME in the raw data

vis$fit$visregFit
sum(vis$fit$visregFit) # still sums to 1 tho

### FIN 