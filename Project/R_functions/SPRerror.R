SPRerror <- function(modpath_vec, Fdyn, niter, mod_names){

	getRE <- function(modpath_vec, iter){
		RE <- rep(NA, length(modpath_vec))
		Rep <- True <- Est <- list()
		for(m in 1:length(modpath_vec)){
			## report file
			Rep[[m]] <- readRDS(file.path(modpath_vec[m], paste0("F_", Fdyn), iter, "Report.rds"))
			
			## estimates of SPR
			Est[[m]] <- with(Rep[[m]], calc_ref(Mat_a=Mat_a, W_a=W_a, M=M, S_a=S_a, R0=exp(beta), F=F_t[length(F_t)], ref=FALSE))
			
			## true values
			if(file.exists(file.path(modpath_vec[m], paste0("F_", Fdyn), iter, "SPR_site.rds"))){
				True_file <- readRDS(file.path(modpath_vec[m], paste0("F_", Fdyn), iter, "SPR_site.rds"))
				True[[m]] <- mean(True_file)
			}
			if(file.exists(file.path(modpath_vec[m], paste0("F_", Fdyn), iter, "SPR_site.rds"))==FALSE){
				True_file <-  readRDS(file.path(modpath_vec[m], paste0("F_", Fdyn), iter, "True.rds"))
				True[[m]] <- True_file$SPR
			}

			RE[m] <- (Est[[m]] - True[[m]])/True[[m]]
		}
		return(RE)
	}

	RE <- t(sapply(1:niter, function(x) getRE(modpath_vec=modpath_vec, iter=x)))
	colnames(RE) <- mod_names

	boxplot(RE, ylim=c(-1,1))
	mtext(side=3, Fdyn, font=2, line=-1)
	abline(h=0, lty=2, col="red")

	return(RE)
}