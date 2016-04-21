# Spatio-temporal models for ecologists
# Week 3 Lab adjusted for Homework 3
# Temporal models

rm(list=ls())
setwd("C:\\Git_Projects\\2016_Spatio-temporal_models\\Week 3 -- Temporal Models\\Homework")

###########################
### load estimation model
###########################

library(TMB)

compile("dlm_hw3.cpp")
dyn.load(dynlib("dlm_hw3"))

###########################
## functions
###########################

DLM_sim <- function(N = 100, seed, sigma.obs = 1, sigma.proc = 0.5, a = -0.5, y1 = 10){
    set.seed(seed)
    ytrue <- numeric(length = N)
    ytrue[1] <- y1
    log.sigma.proc <- log(sigma.proc)
    proc.error <- rnorm(N, mean = 0, sd = sigma.proc)
    log.sigma.obs <- log(sigma.obs)
    for(i in 2:N) {
        ytrue[i] <- a*ytrue[i-1] + proc.error[i-1]
    }
      yobs <- rnorm(N, mean = ytrue, sd = sigma.obs)
    return(list(yobs = yobs, ytrue = ytrue) )
}

DLM_run <- function(sim, model){
		
		#Run TMB model on simulated data
	        N = length(sim$yobs)	
	        data <- list(model=model, y = sim$yobs)
	        parameters <- list(a = 0.5, b=0.5, log_sigma_proc = -1, log_sigma_obs = -1, u = rep(mean(sim$yobs), N))
	        if(model==0){
	        	map <- list()
	        	map[["b"]] <- factor(NA)
	        	obj <- MakeADFun(data, parameters, random = "u", map, DLL = "dlm_hw3")
	        }
	        if(model==1) obj <- MakeADFun(data, parameters, random = "u", DLL = "dlm_hw3")
	        obj$hessian <- FALSE
	        opt <- do.call("optim", obj)
	        sd <- sdreport(obj)	

	        # extract fixed effects:
	  		fixed <- summary(sd, "fixed")

	if(model==0) Outs <- data.frame("a_est"=fixed["a","Estimate"], "proc_est"=exp(fixed["log_sigma_proc", "Estimate"])^2, "obs_est"=exp(fixed["log_sigma_obs", "Estimate"])^2, stringsAsFactors=FALSE)
	if(model==1) Outs <- summary(sd)[which(rownames(summary(sd)) %in% c("a", "b", "sigma_proc", "sigma_obs")),]
	return(Outs)
}

gompertz_sim <- function(N = 100, seed, sigma.obs, sigma.proc, a, b, y1){
	set.seed(seed)
    ytrue <- numeric(length = N)
    ytrue[1] <- y1
    log.sigma.proc <- log(sigma.proc)
    proc.error <- rnorm(N, mean = 0, sd = sigma.proc)
    log.sigma.obs <- log(sigma.obs)
    for(i in 2:N) {
        ytrue[i] <- a + b*ytrue[i-1] + proc.error[i-1]
    }
      yobs <- rnorm(N, mean = ytrue, sd = sigma.obs)
    return(list(yobs = yobs, ytrue = ytrue) )
}



###########################
## simulation experiment
###########################
ar <- c(-0.5, 0, 0.5)
proc <- 0.4
obs <- c(0.5*proc, proc, 2*proc)
iters <- 100

## setupt output
a_output <- list()
a_output[["mean"]] <- a_output[["sd"]] <- a_output[["lower"]] <- a_output[["upper"]] <- matrix(NA, nrow=length(ar), ncol=length(obs))

proc_output <- list()
proc_output[["mean"]] <- proc_output[["sd"]] <- proc_output[["lower"]] <- proc_output[["upper"]] <- matrix(NA, nrow=length(ar), ncol=length(obs))

obs_output <- list()
obs_output[["mean"]] <- obs_output[["sd"]] <- obs_output[["lower"]] <- obs_output[["upper"]] <- matrix(NA, nrow=length(ar), ncol=length(obs))

start <- Sys.time()
## loop over variable AR and observation error values
for(aa in 1:length(ar)){
	for(oo in 1:length(obs)){

		sim_list <- lapply(1:iters, function(x) DLM_sim(a=ar[aa], seed=x, N=100, sigma.obs=obs[oo], sigma.proc=0.4, y1=3))
		res_list <- sapply(1:iters, function(x) DLM_run(sim=sim_list[[x]], model=0))

		a_output[["mean"]][aa,oo] <- mean(unlist(res_list["a_est",]))
		a_output[["sd"]][aa,oo] <- sd(unlist(res_list["a_est",]))
		a_output[["lower"]][aa,oo] <- quantile(unlist(res_list["a_est",]), 0.025)
		a_output[["upper"]][aa,oo] <- quantile(unlist(res_list["a_est",]), 0.975)

		proc_output[["mean"]][aa,oo] <- mean(unlist(res_list["proc_est",]))
		proc_output[["sd"]][aa,oo] <- sd(unlist(res_list["proc_est",]))
		proc_output[["lower"]][aa,oo] <- quantile(unlist(res_list["proc_est",]), 0.025)
		proc_output[["upper"]][aa,oo] <- quantile(unlist(res_list["proc_est",]), 0.975)

		obs_output[["mean"]][aa,oo] <- mean(unlist(res_list["obs_est",]))
		obs_output[["sd"]][aa,oo] <- sd(unlist(res_list["obs_est",]))
		obs_output[["lower"]][aa,oo] <- quantile(unlist(res_list["obs_est",]), 0.025)
		obs_output[["upper"]][aa,oo] <- quantile(unlist(res_list["obs_est",]), 0.975)

		rm(sim_list)
		rm(res_list)
	}
}
end <- Sys.time() - start

saveRDS(a_output, "alpha_output.rds")
saveRDS(proc_output, "process_output.rds")
saveRDS(obs_output, "obs_output.rds")

par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
plot(x=1:3, y=a_output$mean[1,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=a_output$lower[1,],y1=a_output$upper[1,], lwd=2)
abline(h=ar[1], col="red")
mtext(expression("True" ~ alpha ~ "= -0.5"), side=2, line=5, cex=1.5)
mtext("Estimate", side=2, line=3)
mtext(expression(alpha), side=3, line=2, font=2, cex=1.2)

plot(x=1:3, y=proc_output$mean[1,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=proc_output$lower[1,],y1=proc_output$upper[1,], lwd=2)
abline(h=0.4, col="red")
mtext(expression(sigma[p]^2), side=3, line=1, font=2, cex=1.2)

plot(x=1:3, y=obs_output$mean[1,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=obs_output$lower[1,],y1=obs_output$upper[1,], lwd=2)
abline(h=obs[1], col="red")
mtext(expression(sigma[y]^2), side=3, line=1, font=2, cex=1.2)

plot(x=1:3, y=a_output$mean[2,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=a_output$lower[2,],y1=a_output$upper[2,], lwd=2)
abline(h=ar[2], col="red")
mtext(expression("True" ~ alpha ~ "= 0"), side=2, line=5, cex=1.5)
mtext("Estimate", side=2, line=3)

plot(x=1:3, y=proc_output$mean[2,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=proc_output$lower[2,],y1=proc_output$upper[2,], lwd=2)
abline(h=0.4, col="red")

plot(x=1:3, y=obs_output$mean[2,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=obs_output$lower[2,],y1=obs_output$upper[2,], lwd=2)
abline(h=obs[2], col="red")

plot(x=1:3, y=a_output$mean[3,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=a_output$lower[3,],y1=a_output$upper[3,], lwd=2)
axis(1, at=1:3, labels=c(0.2, 0.4, 0.8))
abline(h=ar[3], col="red")
mtext(expression("True" ~ alpha ~ "= 0.5"), side=2, line=5, cex=1.5)
mtext("Estimate", side=2, line=3)

plot(x=1:3, y=proc_output$mean[3,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=proc_output$lower[3,],y1=proc_output$upper[3,], lwd=2)
axis(1, at=1:3, labels=c(0.2, 0.4, 0.8))
abline(h=0.4, col="red")

plot(x=1:3, y=obs_output$mean[3,], pch=16, cex=1.5, xlim=c(0,4), ylim=c(-1,1), xaxt="n", yaxt="n", las=2, ylab="")
segments(x0=1:3,x1=1:3,y0=obs_output$lower[3,],y1=obs_output$upper[3,], lwd=2)
axis(1, at=1:3, labels=c(0.2, 0.4, 0.8))
abline(h=obs[3], col="red")

mtext(side=1, "True Observation Error", line=3, outer=TRUE, cex=1.5)
mtext(side=3, "Parameter", line=4, outer=TRUE, cex=1.5)


###########################
## gompertz model
###########################

gsim <- gompertz_sim(N = 1000, seed=123, sigma.obs=0.2, sigma.proc=0.5, a=1.2, b=0.7, y1=0.5)

par(mfrow=c(2,1), mar=c(1,1,1,1))
plot(x=1, y=1, type="n", ylim=c(0, max(gsim$ytrue)), xlim=c(0,1000), xaxs="i", yaxs="i")
polygon(x=c(0,900,900,0), y=c(0,0,7,7), col="gray")
lines(gsim$ytrue, type="o", lwd=2, ylim = c(min(gsim$ytrue), max(gsim$ytrue)))
points(gsim$yobs, col="blue", pch=19)
lines(gsim$yobs, col="blue", lty=2)

gsim_cut <- lapply(1:length(gsim), function(x) gsim[[x]][901:1000])
names(gsim_cut) <- c("yobs", "ytrue")
plot(gsim_cut$ytrue, type="o", lwd=2, ylim=c(0, max(gsim$ytrue)))
points(gsim_cut$yobs, col="blue", pch=19)
lines(gsim_cut$yobs, col="blue", lty=2)
legend("bottomleft", pch=c(1,19), col=c("black", "blue"), legend=c("state", "observation"), lwd=c(2,1), lty=c(1,2))
mtext(side=1, "Observation #", line=3)


gres <- DLM_run(sim=gsim_cut, model=1)
