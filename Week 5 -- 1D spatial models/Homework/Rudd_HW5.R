rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 5 -- 1D spatial models\\Homework")
require(RandomFields)
require(TMB)

#########################
## simulation function
#########################

Sim_Fn = function( n_i=1000, Scale=2, Sigma2=1, logSD_spatial=0.1, L0=10, Linf_0=100, beta_y=0.02, growth_rate=0.1, mortality_rate=growth_rate*1.6, logSD_resid=0.05){

  # Sample locations
  # n_i = 1000
  y_i = runif( n=n_i, min=32, max=49 )

  # Simulate spatial process
  RMmodel = RMgauss(var=logSD_spatial^2, scale=Scale)
  Linf_i = Linf_0 * exp( RFsimulate(model=RMmodel, x=rep(0,n_i), y=y_i)@data[,1] - Sigma2/2 ) * exp( beta_y*(y_i-40.5) )
  plot( y=Linf_i, x=y_i )

  # Simulate ages of samples
  Survival_a = exp( -mortality_rate * 1:100 )
  a_i = sample( size=n_i, x=1:100, prob=Survival_a, replace=TRUE)
  l_i = Linf_i - (Linf_i-L0) * exp( -growth_rate * a_i ) * exp( rnorm(n_i, mean=-logSD_resid^2/2, sd=logSD_resid) )
  plot( x=a_i, y=l_i )

  # Bundle and return stuff
  DF = data.frame( "y_i"=y_i, "a_i"=a_i, "l_i"=l_i, "Linf_i"=Linf_i)
  return(DF)
}


#############################
## estimation model
#############################

compile( "1d_hw.cpp" )
dyn.load( dynlib("1d_hw") )

#############################
## simulation experiment
#############################

nobs <- 1000
niter <- 100

## simulate data
sim <- Sim_Fn(n_i=nobs)

## data list
Data <- list("Option"=2, "a_i"=sim$a_i, "l_i"=sim$l_i, "loc_i"=sim$y_i)

## parameter inputs
Parameters = list( "beta0"=1, "betaS"=0, "l0"=1, "vbk"=0.2, "ln_sigma2"=log(1), "logit_rho"=1, "epsilon_i"=rep(0, nobs) )

### spatial component
rmse_spatial <- rep(NA, niter)
for(i in 1:niter){
	## run model
	Obj = MakeADFun( data=Data, parameters=Parameters, random="epsilon_i", Map=NULL)	

	Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
	Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
	Report <- Obj$report()

	rmse_spatial[i] <- sqrt((1/nobs)*sum((sim$Linf_i - Report$linf_i)^2))
	
	rm(Report)
	rm(Opt)
	rm(Obj)
}


## no spatial component
### spatial component
rmse_nospatial <- rep(NA, niter)
for(i in 1:niter){

	Map <- list()
	Map[["betaS"]] <- factor(NA)
	## run model
	Obj = MakeADFun( data=Data, parameters=Parameters, random="epsilon_i", Map=Map)	

	Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
	Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
	Report <- Obj$report()

	rmse_nospatial[i] <- sqrt((1/nobs)*sum((sim$Linf_i - Report$linf_i)^2))
	
	rm(Report)
	rm(Opt)
	rm(Obj)
}


Map <- list()
Map[["betaS"]] <- factor(NA)


Sdreport <- sdreport(Obj)

rmse <- sqrt((1/nobs)*sum((sim$Linf_i - Report$linf_i)^2))