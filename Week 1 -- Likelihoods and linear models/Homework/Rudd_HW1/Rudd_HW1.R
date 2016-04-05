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

runModel <- function(data, parameters, map){

	dyn.load( dynlib("glm_hw") )

	Obj <- MakeADFun( data=data, parameters=parameters, map=map)
	Opt <- tryCatch(nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr ), error=function(e) NA)
	if(all(is.na(Opt))==FALSE) Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
	Report <- Obj$report()

	Outs <- NULL
	Outs$Opt <- Opt
	Outs$Report <- Report
	return(Outs)
}

crossValidate <- function(K, data, parameters, map){

   # K <- 10
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

   out <- mean(PredNLL_K / table(Partition_i))

   return(out)
}

simulateData <- function(observations=12210, like, report){
	
	if(like==1) sim <- (1-rbinom(observations, size=1, prob=report$zero_prob))*rlnorm(observations, report$logpred_i, report$sd)

	if(like==2) sim <- (1-rbinom(observations, size=1, prob=report$zero_prob))*rgamma(observations, shape=report$logpred_i, scale=report$sd)

	if(like==3) sim <- rnbinom(observations, mu=report$logpred_i, size=(report$logpred_i^2)/(report$sd^2 - report$logpred_i))

	return(sim)
}


##########################################
### Part 1 - case study demonstration
##########################################
## data
year <- EBS_pollock_data$year
catch <- EBS_pollock_data$catch
# temp <- EBS_pollock_data$waterTmpC
# 	temp[which(temp==-9999)] <- NA
# mean_temp <- mean(temp, na.rm=TRUE)
# temp_dist <- dlnorm(temp, mean_temp, sd(temp, na.rm=TRUE))
# temp_cov <- temp_dist/mean(temp_dist, na.rm=TRUE)
# temp_cov[which(is.na(temp_cov))] <- 1


## design matrix
X <- NULL
X[[1]] <- as.matrix(rep(1, length(catch)))
# X[[2]] <- as.matrix(temp_cov)

## data and parameter lists
Data <- NULL
Data[[1]] <- list("like"=1, "y_i"=catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(catch))) # delta lognormal
Data[[2]] <- list("like"=2, "y_i"=catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(catch))) # delta gamma
Data[[3]] <- list("like"=3, "y_i"=catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(catch))) # negative binomial
# Data[[4]] <- list("like"=1, "y_i"=catch, "X_ij"=as.matrix(X[[2]]), "include"=rep(0,length(catch))) # delta lognormal with temperature covariate


Parameters <- NULL
Parameters[[1]] <- list("b_j"=rep(0, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
Parameters[[2]] <- list("b_j"=rep(1e-5, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
Parameters[[3]] <- list("b_j"=rep(1e-5, ncol(X[[1]])), "theta_z"=c(log(1), log(1)))
# Parameters[[4]] <- list("b_j"=rep(0, ncol(X[[2]])), "theta_z"=c(log(1), log(1)))


Map <- list()
Map[["theta_z"]] <- 1:length(Parameters[["theta_z"]])
Map[["theta_z"]][1] <- NA
Map[["theta_z"]] <- factor(Map[["theta_z"]])


## model 1 - delta lognormal
model1 <- runModel(data=Data[[1]], parameters=Parameters[[1]], map=NULL)
Report1 <- model1$Report
Like1 <- model1$Opt$objective
npar1 <- nrow(model1$Opt$diagnostic)
pprob1 <- crossValidate(K=10, data=Data[[1]], parameters=Parameters[[1]], map=NULL)

## model 2 - delta gamma
model2 <- runModel(data=Data[[2]], parameters=Parameters[[2]], map=NULL)
Report2 <- model2$Report
Like2 <- model2$Opt$objective
npar2 <- nrow(model2$Opt$diagnostic)
pprob2 <- crossValidate(K=10, data=Data[[2]], parameters=Parameters[[2]], map=NULL)


## model 3 - negative binomial
model3 <- runModel(data=Data[[3]], parameters=Parameters[[3]], map=Map)
Report3 <- model3$Report
Like3 <- model3$Opt$objective
npar3 <- nrow(model3$Opt$diagnostic)
pprob3 <- crossValidate(K=10, data=Data[[3]], parameters=Parameters[[3]], map=Map)



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
Sim_catch3 <- simulateData(like=3, report=Report3)
hist(log(1+Sim_catch3), freq=FALSE, add=TRUE, col=rgb(1,1,0,0.2))

# hist(log(1+catch), freq=FALSE, col=rgb(1,0,0,0.2))
# set.seed(123)
# Sim_catch4 <- (1-rbinom(1e5, size=1, prob=Report4$zero_prob))*rlnorm(1e5, Report4$logpred_i, Report4$sd)
# hist(log(1+Sim_catch4), freq=FALSE, add=TRUE, col=rgb(1,1,0,0.2))

##########################################
### Part 2 - simulation experiment
##########################################

like_vec <- 1:3
sim_vec <- 1:100

## iterations by simulation model by estimation model
results <- relerr <- array(NA, dim=c(length(sim_vec),length(like_vec),length(like_vec)))

start <- Sys.time()
for(om in 1:length(like_vec)){
	set.seed(1234)
	for(ss in 1:length(sim_vec)){
		if(like_vec[om]==1) Report <- Report1
		if(like_vec[om]==2) Report <- Report2
		if(like_vec[om]==3) Report <- Report3
		sim_catch <- simulateData(observations=12210, like=like_vec[om], report=Report)


		for(em in 1:length(like_vec)){
			input_data <- list("like"=like_vec[em], "y_i"=sim_catch, "X_ij"=as.matrix(X[[1]]), "include"=rep(0,length(sim_catch)))

			if(like_vec[em]==1){
				input_params <- Parameters[[1]]
				input_map <- NULL
			}
			if(like_vec[em]==2){
				input_params <- Parameters[[2]]
				input_map <- NULL
			}
			if(like_vec[em]==3){
				input_params <- Parameters[[3]]
				input_map <- Map
			}

			model <- runModel(data=input_data, parameters=input_params, map=input_map)
			results[sim_vec[ss], like_vec[om], like_vec[em]] <- model$Report$b_j
			relerr[sim_vec[ss], 1, like_vec[em]] <- (model$Report$b_j - model1$Report$b_j)/model1$Report$b_j
			relerr[sim_vec[ss], 2, like_vec[em]] <- (model$Report$b_j - model2$Report$b_j)/model2$Report$b_j
			relerr[sim_vec[ss], 3, like_vec[em]] <- (model$Report$b_j - model3$Report$b_j)/model3$Report$b_j

			rm(model)
		}

	}
}
end <- Sys.time() - start
