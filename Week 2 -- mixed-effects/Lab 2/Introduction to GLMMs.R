

setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 2 -- mixed-effects\\Lab 2")
Use_REML = FALSE

############
# Generalized linear mixed model
############
library(lme4)

###### Simulate data
# Parameters
Nsite = 10
Nobs_per_site = 10
Site_logMean = log(10)
Site_logSd = 1

# Bookkeeping
s_i = rep( 1:Nsite, each=Nobs_per_site)

# Simulation
z_s = rnorm( Nsite, mean=0, sd=Site_logSd )
Mean_s = exp( Site_logMean + z_s )
y_i = rpois( Nsite*Nobs_per_site, lambda=Mean_s[s_i] )

# Plot data
library(lattice)
histogram( ~ y_i | factor(s_i), breaks=seq( min(y_i), max(y_i), length=10), type="density", panel=function(x,...){ panel.histogram(x, ...); panel.mathdensity(dmath=dnorm, col="black", args = list(mean=mean(x),sd=sd(x))) } )      #

###### Fit using R
	## zero means no intercept
	## factor means you want a value for each size
	## glm automatically assumes canonical link
		## Poisson is easiest to estimate using log link
		## help for GLM will tell you what the canonical links are for each family

# No site level (Not recommended)
GLM = glm( y_i ~ 1, family="poisson" )
print( summary(GLM) )


# Using fixed effects (Not recommended)
	## global intercept
GLM = glm( y_i ~ 0 + factor(s_i), family="poisson" )
print( summary(GLM) )

# Using mixed effects (Recommended) -- doesn't appear to use REML
	## easiest to think about random effects if they are zero-centered, requires fixed effect
library(lme4)
GLMM = glmer( y_i ~ 1 + (1 | factor(s_i)), family="poisson" )
print( summary(GLMM) )

####################
# Fit using TMB
####################
library(TMB)

# Compile model
Version = "glmm"
compile( paste0(Version,".cpp") )

# Build inputs
Data = list( "n_y"=length(y_i), "n_s"=length(unique(s_i)), "s_i"=s_i-1, "y_i"=y_i)
Parameters = list( "x0"=-10, "log_sdz"=2, "z_s"=rep(0,Data$n_s) )
Random = c("z_s")
if( Use_REML==TRUE ) Random = union( Random, "x0")

# Build object
dyn.load( dynlib(Version) )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)  #

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
	## faster when everything is normally distributed
	## newton method works perfectly when you have a perfect quadratic function
	## In this case - likelihood function is skewed normal - overshooting every time
	## takes longer, needs several iterations to optimize random effects given fixed effects
	## actually quasi Newton algorithm - ustep keeps track of curvature relatively to a quadratic likelihood -- ustep telling it not to make as big of steps as it wants to, will go to 1 in a well-behaved model
	## if ustep is not getting to 1, then final gradients are probably not getting low
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  Opt[["total_time"]] = Sys.time() - start_time

# Get reporting and SEs
Report = Obj$report()
  SD = sdreport( Obj ) ## delta method to figure out standard errors (inverse hessian)

  ## standard empirical bayes - just trying to estimate parameters, not random effect
  ## want to make a statement about epsilon - predict them to have a value that maximizes the likelihood function

######################
# Compare estimates
######################

# Global mean
c( fixef(GLMM), Report$x0 )

# Random effects
cbind( "True"=z_s, ranef(GLMM)[['factor(s_i)']], Report$z_s )


