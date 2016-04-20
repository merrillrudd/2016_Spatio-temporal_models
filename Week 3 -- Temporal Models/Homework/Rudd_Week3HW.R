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

ar <- c(-0.5, 0, 0.5)
proc <- 0.4
obs <- c(0.5*proc, proc, 2*proc)
iters <- 100

###########################
## functions
###########################

DLM_sim <- function(N = 100, seed, sigma.obs = 1, sigma.proc = 0.5, a = -0.5, y1 = 10){
    set.seed(seed)
    N <- 100
    ytrue <- numeric(length = N)
    ytrue[1] <- y1
    log.sigma.proc <- log(sigma.proc)
    proc.error <- rnorm(N, mean = 0, sd = sigma.proc)
    log.sigma.obs <- log(sigma.obs)
    for(i in 2:N) {
        ytrue[i] <- a*ytrue[i-1] + proc.error[i-1]
    }
      x <- seq_len(N)
      yobs <- rnorm(N, mean = ytrue, sd = sigma.obs)
    return(list(yobs = yobs, ytrue = ytrue) )
}

DLM_run <- function(sim){
		
		#Run TMB model on simulated data
	        N = length(sim$yobs)	
	        data <- list(y = sim$yobs)
	        parameters <- list(a = 0.5, log_sigma_proc = -1, log_sigma_obs = -1, u = rep(mean(sim$yobs), N))
	        obj <- MakeADFun(data, parameters, random = "u", DLL = "dlm_hw3")
	        obj$hessian <- FALSE
	        opt <- do.call("optim", obj)
	        sd <- sdreport(obj)	

	        # extract fixed effects:
	  		fixed <- summary(sd, "fixed")

	Outs <- data.frame("a_est"=fixed["a","Estimate"], "proc_est"=exp(fixed["log_sigma_proc", "Estimate"])^2, "obs_est"=exp(fixed["log_sigma_obs", "Estimate"])^2, stringsAsFactors=FALSE)
	return(Outs)
}



###########################
## simulation
###########################

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
		res_list <- sapply(1:iters, function(x) DLM_run(sim=sim_list[[x]]))

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




plot(sim$ytrue, type = 'b', col = 2, pch = 's', ylim = c(min(sim$yobs), max(sim$yobs)))
lines(sim$yobs)
points(sim$yobs, pch = 'o')

  # extract estimated process:
  u <- summary(sd, "random")[, "Estimate"]
  u_se <- summary(sd, "random")[, "Std. Error"]

  plot(u, type="l", lwd=2)
  polygon(x=c(1:length(u_se), length(u_se):1), y=c(u-u_se, rev(u)+rev(u_se)), col=rgb(1,0,0,0.5))
  lines(sim$ytrue, lty=2, lwd=3)
  
 
  


  
