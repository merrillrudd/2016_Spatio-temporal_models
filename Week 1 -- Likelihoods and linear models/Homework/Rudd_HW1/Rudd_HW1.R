rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 1 -- Likelihoods and linear models\\Homework\\Rudd_HW1")

# devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
library( SpatialDeltaGLMM )
data( EBS_pollock_data )
library(TMB)

## compile template file
compile( "glm_hw.cpp" )

## data
catch <- EBS_pollock_data$catch

## design matrix
X_list <- NULL
X_list[[1]] <- as.matrix(rep(1, length(catch)))

dyn.load( dynlib("glm_hw") )

## model 1 - delta lognormal
Data <- list(model=1, "y_i"=catch, "X_ij"=as.matrix(X_list[[1]]))
Parameters <- list("b_j"=rep(0,ncol(X_list[[1]])), "theta_z"=c(0,0))
Obj <- MakeADFun( data=Data, parameters=Parameters)

Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics <- data.frame( "name"=names(Obj$par), "Est"=Opt$par, 
	"final_gradient"=as.vector(Obj$gr(Opt$par)) )

SD1 <- sdreport( Obj )

Report1 = Obj$report()

## model 2 - delta gamma 
Data <- list(model=1, "y_i"=catch, "X_ij"=as.matrix(X_list[[1]]))
Parameters <- list("b_j"=rep(0,ncol(X_list[[1]])), "theta_z"=c(0,0))
Obj <- MakeADFun( data=Data, parameters=Parameters)

Opt <- nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr )
Opt$diagnostics <- data.frame( "name"=names(Obj$par), "Est"=Opt$par, 
	"final_gradient"=as.vector(Obj$gr(Opt$par)) )

SD2 <- sdreport( Obj )

Report2 = Obj$report()
