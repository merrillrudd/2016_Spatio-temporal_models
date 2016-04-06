rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 1 -- Likelihoods and linear models\\Homework\\Rudd_HW1")

# devtools::install_github("nwfsc-assess/geostatistical_delta-GLMM")
library( SpatialDeltaGLMM )
data( EBS_pollock_data )
library(TMB)

## compile template file
compile( "glm_hw.cpp" )

##########################################
### functions
##########################################

runModel <- function(data, parameters, map=NULL){

	dyn.load( dynlib("glm_hw") )

	Obj <- MakeADFun( data=data, parameters=parameters, map=map)
	Opt <- tryCatch(nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr), error=function(e) NA)
	if(all(is.na(Opt))==FALSE){
		Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
		
		ii <- 0
	    while(abs(max(Opt$diagnostics$"final_gradient"))>0.05){
    		print("re-running, final gradient too high")
    		ParList <- Obj$env$parList( Obj$env$last.par.best )
                Obj <- MakeADFun(data=data, parameters=ParList) 
                Opt <- tryCatch( nlminb( start= Obj$env$last.par.best*runif(length(Obj$par),0.9,1.1), objective=Obj$fn, gradient=Obj$gr), error=function(e) NA)
                if(all(is.na(Opt))==FALSE) Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
            ii <- ii + 1
            if(ii==25) break
        }
    }
	Report <- Obj$report()

	Outs <- NULL
	Outs$Opt <- Opt
	Outs$Report <- Report
	return(Outs)
}

crossValidate <- function(K, data, parameters, map=NULL){

   Partition_i <- sample(x=1:K, size=length(catch), replace=TRUE)
   PredNLL_K <- rep(NA, K)

   for(k in 1:K){
   	kData <- data[which(names(data) != "include")]
   	kData$include <- ifelse(Partition_i==k,1,0)
   	model <- runModel(data=kData, parameters=parameters, map=map)
   	PredNLL_K[k] <- model$Report$pred_jnll
   	rm(model)
   	rm(kData)
   }

   out <- mean(PredNLL_K / table(Partition_i), na.rm=TRUE)

   return(out)
}

simulateData <- function(observations=12210, like, report){
	
	if(like==1) sim <- (1-rbinom(observations, size=1, prob=report$zero_prob))*rlnorm(observations, report$logpred_i, report$sd)

	if(like==2) sim <- (1-rbinom(observations, size=1, prob=report$zero_prob))*rgamma(observations, shape=1/(report$sd^2), scale=exp(report$logpred_i)*(report$sd^2))

	return(sim)
}


##########################################
### Part 1 - case study demonstration
##########################################

## data
year <- EBS_pollock_data$year
catch <- EBS_pollock_data$catch


## design matrix
X <- NULL
X[[1]] <- as.matrix(rep(1, length(catch)))
X[[2]] <- as.matrix(cbind(rep(1,length(catch)), EBS_pollock_data$lat, EBS_pollock_data$long))

## data and parameter lists
Data <- NULL
Data[[1]] <- list("like"=1, "y_i"=catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(catch))) # delta lognormal
Data[[2]] <- list("like"=2, "y_i"=catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(catch))) # delta gamma
Data[[3]] <- list("like"=2, "y_i"=catch, "X_ij"=as.matrix(X[[2]]), "include"=rep(0,length(catch))) # delta lognormal with lat/lon covariate


Parameters <- NULL
Parameters[[1]] <- list("b_j"=rep(0, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
Parameters[[2]] <- list("b_j"=rep(0, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
Parameters[[3]] <- list("b_j"=rep(0, ncol(X[[2]])), "theta_z"=c(log(1), log(1)))


## model 1 - delta lognormal
model1 <- runModel(data=Data[[1]], parameters=Parameters[[1]])
Report1 <- model1$Report
Opt1 <- model1$Opt
Like1 <- model1$Opt$objective
npar1 <- nrow(model1$Opt$diagnostic)
set.seed(1234)
pprob1 <- crossValidate(K=10, data=Data[[1]], parameters=Parameters[[1]])

## model 2 - delta gamma
model2 <- runModel(data=Data[[2]], parameters=Parameters[[2]])
Report2 <- model2$Report
Opt2 <- model2$Opt
Like2 <- model2$Opt$objective
npar2 <- nrow(model2$Opt$diagnostic)
set.seed(1234)
pprob2 <- crossValidate(K=10, data=Data[[2]], parameters=Parameters[[2]])


## model 3 - negative binomial
model3 <- runModel(data=Data[[3]], parameters=Parameters[[3]])
Report3 <- model3$Report
Opt3 <- model3$Opt
Like3 <- model3$Opt$objective
npar3 <- nrow(model3$Opt$diagnostic)
set.seed(1234)
pprob3 <- crossValidate(K=10, data=Data[[3]], parameters=Parameters[[3]])



## plot model fits
par(mfrow=c(2,2))
hist(log(1+catch), freq=FALSE, col=rgb(1,0,0,0.2))
set.seed(123)
Sim_catch1 <- simulateData(like=1, report=Report1)
hist(log(1+Sim_catch1), freq=FALSE, add=TRUE, col=rgb(0,0,1,0.2))

hist(log(1+catch), freq=FALSE, col=rgb(1,0,0,0.2))
set.seed(123)
Sim_catch2 <- simulateData(like=2, report=Report2)
hist(log(1+Sim_catch2), freq=FALSE, add=TRUE, col=rgb(0,1,0,0.2))

hist(log(1+catch), freq=FALSE, col=rgb(1,0,0,0.2))
set.seed(123)
Sim_catch3 <- simulateData(like=2, report=Report3)
hist(log(1+Sim_catch3), freq=FALSE, add=TRUE, col=rgb(1,1,0,0.2))

##########################################
### Part 2 - simulation experiment
##########################################

sim_vec <- 1:100

estimates <- relerr <- array(NA, dim=c(length(sim_vec), 3, 3))

start <- Sys.time()
for(ss in sim_vec){

	for(om in 1:3){
		set.seed(ss)
		if(om==1) simdata <- simulateData(like=1, report=Report1) ## data truly arise from delta lognormal process 
		if(om==2) simdata <- simulateData(like=2, report=Report2) ## data truly arise from delta gamma process
		if(om==3) simdata <- simulateData(like=2, report=Report3) ## data truly arise from delta gamma process with covariates

	    for(em in 1:3){
    		## assuming data arise from lognormal process
    		if(em==1){
    			inputdata <- list("like"=1, "y_i"=simdata, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(simdata)))
    			inputparams <- list("b_j"=rep(0, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
    		}
    		## assuming data arise from gamma process
    		if(em==2){
    			inputdata <- list("like"=2, "y_i"=simdata, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(simdata)))
    			inputparams <- list("b_j"=rep(0, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
    		}
    		## assuming data arise from gamma process with covariates
    		if(em==3){
    			inputdata <- list("like"=2, "y_i"=simdata, "X_ij"=as.matrix(X[[2]]), "include"=rep(0,length(simdata)))
    			inputparams <- list("b_j"=rep(0, ncol(X[[2]])), "theta_z"=c(log(1), log(1)))
    		}
    		run <- runModel(data=inputdata, parameters=inputparams)
    		estimates[ss,em,om] <- run$Report$b_j[1]
    		if(om==1) relerr[ss,em,om] <- (run$Report$b_j[1] - Report1$b_j[1])/Report1$b_j[1]
    		if(om==2) relerr[ss,em,om] <- (run$Report$b_j[1] - Report2$b_j[1])/Report2$b_j[1]
    		if(om==3) relerr[ss,em,om] <- (run$Report$b_j[1] - Report3$b_j[1])/Report3$b_j[1]
    	}
    }
}
end <- Sys.time() - start
