rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 2 -- mixed-effects\\Homework\\Rudd_HW2")
library(TMB)

######################
## compile exe
######################

compile("glmm_hw.cpp")

######################
## settings
######################

nsites <- 10	## number of sites
nobs <- 10 		## number of observations per site
tobs <- nsites*nobs

mu <- 2 		## log-mean for expected counts
var_s <- 1 		## variance of among-site variability
var_y <- 0.5	## variance of overdispersion
sig_s <- sqrt(var_s) ## standard deviation - among site variability
sig_y <- sqrt(var_y) ## standard deviation - overdispersion

#############################
## data generation
#############################
set.seed(1)
ldens <- rnorm(nsites, mu, sig_s)

set.seed(2)
lobs <- unlist(lapply(1:nsites, function(x) rnorm(nobs, ldens[x], sig_y)))

set.seed(3)
sim_counts <- rpois(tobs, exp(lobs))
sim_index <- as.vector(sapply(1:nsites, function(x) rep(x-1,nobs)))

##############################
## setup estimation models
##############################

dyn.load( dynlib("glmm_hw") )

## setup input lists
Data <- lapply(1:4, function(x) list(model=x, tobs=tobs, nsites=nsites, y_j=sim_counts, site_j=sim_index))

Parameters <- list("mu"=mu, "logsig_s"=log(sig_s), "logsig_y"=log(sig_y), "tau"=rep(0,nsites), "epsilon"=rep(0, tobs))

Random <- list()
Random[[1]] <- NULL
Random[[2]] <- c("tau")
Random[[3]] <- c("epsilon")
Random[[4]] <- c("tau", "epsilon")

Map <- list()
Map[[1]] <- list("logsig_s"=factor(NA), "logsig_y"=factor(NA), "tau"=factor(rep(NA, length(Parameters$tau))), "epsilon"=factor(rep(NA, length(Parameters$epsilon))))
Map[[2]] <- list("logsig_y"=factor(NA), "epsilon"=factor(rep(NA, length(Parameters$epsilon))))
Map[[3]] <- list("logsig_s"=factor(NA), "tau"=factor(rep(NA, length(Parameters$tau))))

runModel <- function(data, parameters, random, map){
	Obj <- MakeADFun(data=data, parameters=parameters, random=random, map=map)
	Opt <- nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr)
	Opt$diagnostics <- data.frame("name"=names(Obj$par), "Est"=Opt$par, "final_gradient"=as.vector(Obj$gr(Opt$par)))
	Report <- Obj$report()
	Sdreport <- sdreport(Obj)

	Outs <- NULL
	Outs$Opt <- Opt
	Outs$Report <- Report
	Outs$Sdreport <- Sdreport
	return(Outs)
}

## run models 
model1 <- runModel(data=Data[[1]], parameters=Parameters, random=Random[[1]], map=Map[[1]])

model2 <- runModel(data=Data[[2]], parameters=Parameters, random=Random[[2]],map=Map[[2]])

model3 <- runModel(data=Data[[3]], parameters=Parameters, random=Random[[3]], map=Map[[3]])

model4 <- runModel(data=Data[[4]], parameters=Parameters, random=Random[[4]],map=NULL)

### plot predicted vs. observed counts
par(mfrow=c(2,2))
plot(sim_index, sim_counts, pch=16)
points(sim_index, model1$Report$pred_j, col="blue", lwd=2, cex=1.2)

plot(sim_index, sim_counts, pch=16)
points(sim_index, model2$Report$pred_j, col="blue", lwd=2, cex=1.2)

plot(sim_index, sim_counts, pch=16)
points(sim_index, model3$Report$pred_j, col="blue", lwd=2, cex=1.2)

plot(sim_index, sim_counts, pch=16)
points(sim_index, model4$Report$pred_j, col="blue", lwd=2, cex=1.2)


#############################
## simulation testing
#############################

## generate 100 datasets
nsim <- 100
set.seed(3)
simdata <- lapply(1:nsim, function(x) rpois(tobs, exp(lobs)))

est_mu <- array(NA, dim=c(4, 2, nsim))
	colnames(est_mu) <- c("Estimate", "SE")
for(i in 1:nsim){

	Data_sim <- lapply(1:4, function(x) list(model=x, tobs=tobs, nsites=nsites, y_j=simdata[[i]], site_j=sim_index))

	## run models 
    model1 <- runModel(data=Data_sim[[1]], parameters=Parameters, random=Random[[1]], map=Map[[1]])    

    model2 <- runModel(data=Data_sim[[2]], parameters=Parameters, random=Random[[2]],map=Map[[2]])    

    model3 <- runModel(data=Data_sim[[3]], parameters=Parameters, random=Random[[3]], map=Map[[3]])    

    model4 <- runModel(data=Data_sim[[4]], parameters=Parameters, random=Random[[4]],map=NULL)

    est_mu[1,,i] <- summary(model1$Sdreport)["mu",]
    est_mu[2,,i] <- summary(model2$Sdreport)["mu",]
    est_mu[3,,i] <- summary(model3$Sdreport)["mu",]
    est_mu[4,,i] <- summary(model4$Sdreport)["mu",]

}

par(mfrow=c(2,2), mar=c(0,0,0,0), omi=c(1,1,1,1))
plot(est_mu[1,1,], ylim=c(1.6, 2.9), xaxt="n", yaxt="n", pch=19, xaxs="i", yaxs="i")
segments(x0=1:nsim, x1=1:nsim, y0=est_mu[1,1,]-1.96*est_mu[1,2,], y1=est_mu[1,1,]+1.96*est_mu[1,2,])
count1 <- length(which(est_mu[1,1,]-1.96*est_mu[1,2,] < mu & est_mu[1,1,]+1.96*est_mu[1,2,] > mu))
abline(h=mu, col="red", lwd=2)
axis(2, pretty(c(1.6, 2.9))[-1])
mtext(side=1, line=-2, "No site or \noverdispersion effects", font=2)

plot(est_mu[2,1,], ylim=c(1.6, 2.9), xaxt="n", yaxt="n", pch=19, xaxs="i", yaxs="i")
segments(x0=1:nsim, x1=1:nsim, y0=est_mu[2,1,]-1.96*est_mu[2,2,], y1=est_mu[2,1,]+1.96*est_mu[2,2,])
count2 <- length(which(est_mu[2,1,]-1.96*est_mu[2,2,] < mu & est_mu[2,1,]+1.96*est_mu[2,2,] > mu))
abline(h=mu, col="red", lwd=2)
mtext(side=1, line=-2, "Among-site variability", font=2)


plot(est_mu[3,1,], ylim=c(1.6, 2.9), pch=19, xaxs="i", yaxs="i")
segments(x0=1:nsim, x1=1:nsim, y0=est_mu[3,1,]-1.96*est_mu[3,2,], y1=est_mu[3,1,]+1.96*est_mu[3,2,])
count3 <- length(which(est_mu[3,1,]-1.96*est_mu[3,2,] < mu & est_mu[3,1,]+1.96*est_mu[3,2,] > mu))
abline(h=mu, col="red", lwd=2)
mtext(side=3, line=-2, "Overdispersion", font=2)


plot(est_mu[4,1,], ylim=c(1.6, 2.9), yaxt="n", pch=19, xaxs="i", yaxs="i")
segments(x0=1:nsim, x1=1:nsim, y0=est_mu[4,1,]-1.96*est_mu[4,2,], y1=est_mu[4,1,]+1.96*est_mu[4,2,])
count4 <- length(which(est_mu[4,1,]-1.96*est_mu[4,2,] < mu & est_mu[4,1,]+1.96*est_mu[4,2,] > mu))
abline(h=mu, col="red", lwd=2)
mtext(side=3, line=-3, "Among-site variability \nand overdispersion", font=2)

mtext(side=1, "Simulated Dataset", outer=TRUE, line=3)
mtext(side=2, "Estimate", outer=TRUE, line=3)
