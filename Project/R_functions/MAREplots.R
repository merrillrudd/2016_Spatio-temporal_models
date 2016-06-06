MAREplots <- function(F_scenarios, R_scenarios, data_scenarios, SigR, 
	itervec, param, rel_err=TRUE, simdir, lh_num, lh_dat, Ncomp_scenarios=FALSE,
  sens_param=FALSE, sens_val=FALSE, RecType, OMselex, EMselex, absolute){

  ncol <- length(F_scenarios)*length(data_scenarios)
  nrow <- length(R_scenarios)

  if(RecType==0) RecTypeDesc <- "shared"
  if(RecType==1) RecTypeDesc <- "BH"
  if(RecType==2) RecTypeDesc <- "constant"

  par(mfrow=c(nrow, ncol), mar=c(0,0,0,0), omi=c(1,1,1,1))

  if(simdir=="base"){
  	for(f in 1:length(F_scenarios)){
  		for(r in 1:length(R_scenarios)){
  			for(d in 1:length(data_scenarios)){
				mare <- are <- mre <- re <- rep(NA, length=length(itervec))  
				x <- file.path(test_dir, "base", paste0("LH_", lh_num), 
          paste0("OM", OMselex, "_EM", EMselex), paste0("RecType_", RecTypeDesc),
          paste0("SigR_", SigR), paste0("Fdyn_", F_scenarios[f]), paste0("Rdyn_", R_scenarios[r]),
          data_scenarios[d])
				for(iter in 1:length(itervec)){
					if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
	            	True <- readRDS(file.path(x, iter, "True.rds"))
				    if(param=="FFref"){
              pred <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                 M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.35)$root
              true <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=True$Mat_a, W_a=True$W_a, 
                 M=True$M, S_a=True$S_a, R0=True$R0, ref=0.35)$root
            }
            if(param=="Depl"){
              pred <- Report$Depl[length(Report$Depl)]
              true <- True$D_t[length(True$D_t)]
            }
              re[iter] <- (pred - true)/true
              mre[iter] <- median(c(re[1:iter]))

              are[iter] <- abs(pred - true)/true
              mare[iter] <- median(c(are[1:iter]))


          }
            if(absolute==FALSE){
                  plot(x=itervec, y=mre, xaxt="n", yaxt="n", pch=19, col="#AAAAAAAA", ylim=c(-1,1))
                  abline(h=0, lty=2)
            }
            if(absolute==TRUE){
                  plot(x=itervec, y=mare, xaxt="n", yaxt="n", pch=19, col="#AA0000AA", ylim=c(-1,1))
                  abline(h=0, lty=2)
            }
  			}
  		}
  	}
  }
}
