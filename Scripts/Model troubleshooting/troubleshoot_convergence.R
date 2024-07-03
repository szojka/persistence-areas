
## When there are convergence issues:

#1.  https://mran.microsoft.com/snapshot/2017-11-04/web/packages/glmmTMB/vignettes/troubleshooting.html
# problem is non-positive hessian matrix. Sometimes because model is overpareratized, like in the zi term.
# zi might be due to a different term that's not in the model. In my case, that could be biotic vs abiotic
# treatment. Try fitting model with A and B in them as interaction w/ x, and zi ~ treatment.
# this works - no convergence issues

#2.  #mp1a <- update(mp1a, start=getME(mp1a, c("theta")), control=glmerControl(optimizer="Nelder_Mead")) # doesn't help

#3.  compois distribution can help


# for non-positive hessian matrix 
# always scale x (in ggpredict common solution to confidence interavals not coming out)
# sometimes need to round to 2 (if using ggpredict) for connecting for backtransform 

#############################################################################
# FOR MCLOGIT MULTINOMIAL:
# 
# # troubleshoot convergence
# # BEN BOLKER: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2013q1/019968.html
# mfit <- mblogit(y_mat ~ green_index_scaled,
#                 random = c(~ 1|site/grid, ~1|block),
#                 weights = replicates, 
#                 data = green_all_DLMESS,
#                 control = mmclogit.control(maxit = 1000,
#                 epsilon = 1e-08, trace=TRUE)) # converged on last iteration

# WHAT I HAVE TRIED TO FIX ERROR INNER FALSE CONVERGENCE (8)
# change method for random effects (MQL doesnt' work and PQL has same issue)
# estimator changed to ML or REML
# avoid.increase = FALSE / TRUE
# break.on.increase = FALSE / TRUE
# break.on.infinite	= TRUE / FALSE
# break.on.negative	= FALSE / TRUE
# trace=FALSE no convergence anywhere
# epsilon = 1e-09, 1e-07, 1e-06, 1e-05, 1e-04, 1e-10, 1e-11
# catCov= 'diagonal', # nlminb message:singular convergence (7)
# 'free', false convergence (8)
# 'single'Warning: Inner iterations did not coverge - nlminb message: singular convergence (7)  https://rdrr.io/cran/mclogit/man/mblogit.html
# separate fitting or A vs B! Remove block. Still false convergence (8)
# simplified random effects: (1|site) no help, (1|grid)no help, (1|block) no help


# this sort of optimizer control doesn't work control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))

#  random slopes do not fix convergence

# Iterations 1 through 4 Warning: Inner iterations did not coverge - nlminb message: false convergence (8) 
#https://stackoverflow.com/questions/40039114/r-nlminb-what-does-false-convergence-actually-mean#:~:text=p.%205%3A%20false%20convergence%3A%20the%20gradient%20%E2%88%87f%20%28x%29,%28XFTOL%29%20%E2%80%94%20V%20%2834%29%20is%20the%20false-convergence%20tolerance.

# Iteration 5 - deviance = 564.7487 - criterion = 1.141349e-11
# converged

# Is it still a problem that the in bw iterations did not convergence?
# Ben Bolker: To eliminate some ‘‘false convergence’’ messages and useless function evaluations, 
# it is necessary to increase the stopping tolerances and, when finite-difference derivative 
# approximations are used, to increase the step-sizes used in estimating derivatives.
# Megan: did not successfully eliminate convergence issues . see maxit comment above.

# Each iteration has nearly identical deviance criterions with or without random effects. Going to say it is appropriate to remove random effects.
