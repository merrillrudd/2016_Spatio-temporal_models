
library( devtools )
library( TMB )
# install_github("nwfsc-assess/geostatistical_delta-GLMM")
library( SpatialDeltaGLMM )

setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 1 MBR\\Lecture 1")

############
# Example 1 -- average CPUE for canary rockfish
############

# Data
data(WCGBTS_Canary_example)
CPUE = WCGBTS_Canary_example$HAUL_WT_KG
# Visualize
png( file="Canary_histogram.png", width=4, height=4, res=200, units="in")
  par( mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
  hist( log(1+CPUE) )
dev.off()
# Probability for example parameters
dnorm( CPUE, mean=20, sd=2 )           # Likelihood for each datum
sum(dnorm( CPUE, mean=20, sd=2, log=TRUE ))     # Log-likelihood for all data

######### Method 1 -- Pre-made functions in R
## estimates 2 parameters - mean and standard deviation
Lm = lm( CPUE ~ 1 )
summary(Lm)


######### Method 2 -- Optimize using R
# Step 1 -- define function
## takes vector of parameters, data
NegLogLike_Fn = function(Par, Data){
  # Parameters
  Mean_hat = Par[1]
  SD_hat = Par[2]
  # Log-likelihood
  LogLike_i = dnorm( Data$y, mean=Mean_hat, sd=SD_hat, log=TRUE )
  NegLogLike = -1 * sum(LogLike_i)
  return( NegLogLike )
}
# step 2 -- minimize negative log-likelihood to estimate parameters
## tagged list - input into TMB
Data = list( 'y'=CPUE )
Start = c(1,1)
NegLogLike_Fn(Par=Start, Data=Data)
Opt = optim( par=Start, fn=NegLogLike_Fn, Data=Data, lower=c(0.01,0.01), upper=Inf, method="L-BFGS-B", hessian=TRUE )
print( Opt$par ) # Estimated parameters
print( sqrt( diag( solve(Opt$hessian) ) ) )# standard errors (HESSIAN used to get asymptotic standard errors)
	## solve command is taking inverse of hessian matrix

###### Method 3 -- Optimize using TMB
# Step 1 -- make and compile template file
compile( "linear_model_v1.cpp" )
  ## creates DLL == dynamically linked library, memory allocation outside of R that efficiently runs programs

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v1") )
  ## TMB:::getUserDLL() <-- tells users what is loaded
Params = list("mean"=0, "log_sd"=0) 
  ## log_sd - likelihood profile of standard deviations and variances is often curved (looks like lognormal distribution), profile of log-variance will be more u-shaped and standard error calculations will be easier to calculate
Data = list( "y_i"=CPUE )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v1")


# Step 3 -- test and optimize
Obj$fn( Obj$par ) ## if we give it possible parameter values it spits out negative log likelihood -- identical to negloglike_fn (same usage), but what's new is we didn't previously have access to the gradients
Obj$gr( Obj$par ) ## get a vector of gradients, where the vector has same length at the number of parameters
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr ) ## can use any nonlinear optimizers in R, pass starting values for parameters, function to calculate negative log likelihood, and function that calculates gradients
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors

############
# Example 1 -- average CPUE for canary rockfish
############

# Define a new design matrix
X = cbind( "CA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD<42,1,0), "OR"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=42&WCGBTS_Canary_example$BEST_LAT_DD<46,1,0), "WA"=ifelse(WCGBTS_Canary_example$BEST_LAT_DD>=46,1,0), "Pass"=WCGBTS_Canary_example$PASS-1.5 )

# Step 1 -- make and compile template file
compile( "linear_model_v2.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("linear_model_v2") )
Params = list("b_j"=rep(0,ncol(X)), "log_sd"=0)
Data = list( "y_i"=CPUE, "X_ij"=X )
Obj = MakeADFun( data=Data, parameters=Params, DLL="linear_model_v2")

# Step 3 -- test and optimize
Obj$fn( Obj$par )
Obj$gr( Obj$par )
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics = data.frame( "name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
Opt$par # estimated parameters
SD = sdreport( Obj ) # standard errors




