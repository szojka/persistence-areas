require(dplyr)
require(glmmTMB)
require(car)
require(visreg)

# So I just wrote some fake code (attached) that roughly mirrors your design – i.e., there are 4 species in each plot, and possibly 4 categories of dispersal limitation, mass effects, etc. If logistic regression is not totally screwing up, for each species, the sum of the fitted probabilities for the 4 categories should be 1 (on average they should also be 0.25 each in this example I’ve made but which category gets a 1 or a 0 is determined by random draw, so there’s some variability there).  
# 
# 
# 
# The fitted probabilities does in fact sum to 1 across the 4 categories, so that’s comforting. Technically the test is asking if the probabilities are different, as opposed to each differing from an expected proportion of 0.25, but if indeed they all had the same probabilities, these probabilities would necessarily be 0.25, so I think conceptually it’s doing the same thing.

#4 species, 10 plots
df<-data.frame(plot=as.factor(rep(1:10,each=16)), species=as.factor(rep(rep(1:4,each=4),10)),type=as.factor(rep(1:4,40)))
df$code<-paste(df$plot,df$species,sep="_")
df$random<-paste(df$plot,df$species,df$type,sep="_")

tmp2<-data.frame()

  
  for(i in 1:length(unique(df$code))) {
    
  tmp<-filter(df,code==unique(code)[i])  
  print(tmp)
  
  tmp1<-sample_n(tmp, 1)
  tmp1$data<-1
  tmp2<-rbind(tmp2,tmp1)
    
  }


df<-left_join(df,tmp2)
df$data<-ifelse(is.na(df$data)==TRUE,0,1)

###########
#stats
###########

lm<-glmmTMB(data ~ type*species + (1|random),data=df,family="binomial")
Anova(lm)
vis<-visreg(lm,xvar="type",scale="response")

vis$fit$visregFit
sum(vis$fit$visregFit) #should sum to about 1
