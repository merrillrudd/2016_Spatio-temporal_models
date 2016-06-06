mean_rel_error <- function(Fscen, Rscen, data_scenarios, SigR, 
	itervec, simdir, lh_num, lh_dat, Ncomp_scenarios=FALSE, param,
  sens_param=FALSE, sens_val=FALSE, RecType, OMselex, EMselex, absolute=TRUE){


    if(RecType==0) RecTypeDesc <- "shared"
    if(RecType==1) RecTypeDesc <- "BH"
    if(RecType==2) RecTypeDesc <- "constant"


	if(simdir=="base"){
			
			mre <- mare <- rep(NA, length(data_scenarios))
      are <- re <- matrix(NA, nrow=length(itervec), ncol=length(data_scenarios))

	    for(dat in 1:length(data_scenarios)){
    		x <- file.path(test_dir, "base", paste0("LH_", lh_num),
    			paste0("OM", OMselex, "_EM", EMselex), 
    			paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR),
    			paste0("Fdyn_", Fscen), paste0("Rdyn_", Rscen),
    			data_scenarios[dat])
    		for(iter in 1:length(itervec)){
       		  if(file.exists(file.path(x, iter, "True.rds"))) True <- readRDS(file.path(x, iter, "True.rds"))
       		  if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
       		  if(file.exists(file.path(x, iter, "True.rds"))==FALSE) next
       		  if(file.exists(file.path(x, iter, "Report.rds"))==FALSE) next

          if(param=="Depl"){
       		  t <- True$D_t[length(True$D_t)]
       		  r <- Report$Depl[length(Report$Depl)]
       		  are[iter,dat] <- abs(r - t)/t
            re[iter,dat] <- (r - t)/t
          }
          if(param=="FFref"){
              rref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                 M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
              tref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                 M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
              rval <- True$F_t[length(True$F_t)]
              tval <- Report$F_t[length(Report$F_t)]
              r <- rval/rref
              t <- tval/tref
              are[iter,dat] <- abs(r - t)/t
              re[iter,dat] <- (r - t)/t
          }
    		}
    		mre[dat] <- median(re[,dat], na.rm=TRUE)
        mare[dat] <- median(are[,dat], na.rm=TRUE)
    	}
    }

    if(simdir=="explore_comp"){

      mre <- mare <- rep(NA, length(Ncomp_scenarios))
      re <- are <- matrix(NA, nrow=length(itervec), ncol=length(Ncomp_scenarios))

	    for(dat in 1:length(Ncomp_scenarios)){
    		x <- file.path(test_dir, "explore_comp", paste0("Nyrs_comp_", Ncomp_scenarios[dat]),
    		    paste0("LH_", lh_num),
    			paste0("OM", OMselex, "_EM", EMselex), 
    			paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR),
    			paste0("Fdyn_", Fscen), paste0("Rdyn_", Rscen),
    			"Poor_Comp_LC")
    		for(iter in 1:length(itervec)){
       		  if(file.exists(file.path(x, iter, "True.rds"))) True <- readRDS(file.path(x, iter, "True.rds"))
       		  if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
       		  if(file.exists(file.path(x, iter, "True.rds"))==FALSE) next
       		  if(file.exists(file.path(x, iter, "Report.rds"))==FALSE) next

          if(param=="Depl"){
       		  t <- True$D_t[length(True$D_t)]
       		  r <- Report$Depl[length(Report$Depl)]
       		  are[iter,dat] <- abs(r - t)/t
            re[iter,dat] <- abs(r - t)/t
          }
          if(param=="FFref"){
              rref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                 M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
              tref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                 M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
              rval <- True$F_t[length(True$F_t)]
              tval <- Report$F_t[length(Report$F_t)]
              r <- rval/rref
              t <- tval/tref
              are[iter,dat] <- abs(r - t)/t
              re[iter,dat] <- (r - t)/t
          }
    		}
    		mre[dat] <- median(re[,dat], na.rm=TRUE)
        mare[dat] <- median(are[,dat], na.rm=TRUE)
    	}
    }

    Outs <- NULL
    Outs$mre <- mre
    Outs$mare <- mare
    Outs$re <- re
    Outs$are <- are
    return(Outs)



}