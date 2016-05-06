rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 5 -- 1D spatial models\\Homework")
require(RandomFields)
require(TMB)

#########################
## simulation function
#########################
######
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

run_model <- function(sim, nobs, spatial){
	## simulate data

	## data list
	Data <- list("Option"=2, "a_i"=sim$a_i, "l_i"=sim$l_i, "loc_i"=sim$y_i)

	## parameter inputs
	Parameters = list( "beta0"=1, "betaS"=0, "l0"=1, "vbk"=0.2, "ln_sigma2"=log(1), "logit_rho"=1, "epsilon_i"=rep(0, nobs) )

	## run model
	if(spatial==TRUE) Obj = MakeADFun( data=Data, parameters=Parameters, random="epsilon_i", map=NULL)	
	if(spatial==FALSE){
		Map <- list()
		Map[["betaS"]] <- factor(NA)

		## run model
		Obj = MakeADFun( data=Data, parameters=Parameters, random="epsilon_i", map=Map)	
	}

	Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
	Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
	Report <- Obj$report()

	rmse <- sqrt((1/nobs)*sum((sim$Linf_i - Report$linf_i)^2))
	return(rmse)
}


#######################
## estimation model
#############################

compile( "1d_hw.cpp" )
dyn.load( dynlib("1d_hw") )

#############################
## simulation experiment
#############################
start <- Sys.time()
nobs <- 1000
niter <- 100
set.seed(123)
rmse_spatial <- rmse_nonspatial <- rep(NA, niter)
for(i in 1:niter){
	sim <- Sim_Fn(n_i=nobs)

	rmse_spatial[i] <- run_model(sim=sim, nobs=nobs, spatial=TRUE)
	rmse_nonspatial[i] <- run_model(sim=sim, nobs=nobs, spatial=FALSE)
}
end <- Sys.time() - start

rmse1_adj <- rmse_spatial[1:20]
rmse2_adj <- rmse_nonspatial[1:20]
boxplot(rmse1_adj, rmse2_adj, ylim=c(-0.1, max(c(rmse1_adj, rmse2_adj))))
abline(h=0, lty=2)
axis(1, at=c(1,2), labels=c("Spatial trend", "No spatial trend"))
mtext("Model", side=1, line=3)
mtext("Root Mean Squared Error", side=2, line=3)