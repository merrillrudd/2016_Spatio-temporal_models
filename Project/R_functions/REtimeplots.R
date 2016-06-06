REtimeplots <- function(F_scenarios, R_scenarios, data_scenarios, 
	SigR, itervec, param, simdir, lh_num, lh_dat, Ncomp_scenarios=FALSE,
	sens_param=FALSE, sens_val=FALSE, RecType, OMselex, EMselex){


  if(simdir=="base"){
    data_vec <- data_scenarios
    data_scenario_names <- rep(NA, length(data_vec))
    data_scenario_names[grep("Rich", data_vec)] <- "Rich"
    data_scenario_names[grep("Moderate_Full", data_vec)] <- "Moderate_Full"
    data_scenario_names[grep("Moderate_Sample", data_vec)] <- "Moderate_Sample"
    data_scenario_names[grep("Poor_Index", data_vec)] <- "Poor_Index"
    data_scenario_names[grep("Poor_Comp", data_vec)] <- "Poor_Comp_Only"
    data_scenario_names[grep("Poor_Rel", data_vec)] <- "Poor_Relative"
  }
  if(simdir=="explore_comp"){
    data_vec <- Ncomp_scenarios
    data_scenario_names <- Ncomp_scenarios
  }

	nrow <- ifelse(length(data_vec)==1, length(R_scenarios), length(data_vec))
	ncol <- ifelse(length(F_scenarios)==1, length(R_scenarios), length(F_scenarios))






	calc_re <- function(F_choose, R_choose, param, simdir, dat_choose){
     if(simdir=="base"){
     	re <- matrix(NA, nrow=length(itervec), ncol=20)

        if(RecType==0) RecTypeDesc <- "shared"
        if(RecType==1) RecTypeDesc <- "BH"
        if(RecType==2) RecTypeDesc <- "constant"

              x <- file.path(test_dir, "base", paste0("LH_", lh_num), paste0("OM", OMselex, "_EM", EMselex), 
                 paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR), paste0("Fdyn_", F_choose), paste0("Rdyn_", R_choose),
                 dat_choose)

          for(iter in 1:length(itervec)){
            if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
            True <- readRDS(file.path(x, iter, "True.rds"))

              if(param=="SPR"){
                pred <- calc_SPR(SB0=Report$SB0, R0=exp(Report$beta), SB_t=Report$SB_t_hat, R_t=Report$R_t_hat)
                true <- calc_SPR(SB0=True$SB0, R0=lh_dat$R0, SB_t=True$SB_t, R_t=True$R_t)                
          	  	re[iter,] <- (pred - true)/true
          	  }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,] <- ((Report$F_t/pred) - (True$F_t/true))/(True$F_t/true)
              } 
              if(param=="Depl"){
                pred <- Report$Depl
                true <- True$D_t
                re[iter,] <- (pred - true)/true
              }
          	
          }
       } ## end base

      if(simdir=="explore_comp"){
      re <- matrix(NA, nrow=length(itervec), ncol=20)

        if(RecType==0) RecTypeDesc <- "shared"
        if(RecType==1) RecTypeDesc <- "BH"
        if(RecType==2) RecTypeDesc <- "constant"

            x <- file.path(test_dir, "explore_comp", paste0("Nyrs_comp_", dat_choose), 
                 paste0("LH_", lh_num),  paste0("OM", OMselex, "_EM", EMselex), paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR), 
                 paste0("Fdyn_", F_choose), paste0("Rdyn_", R_choose),
                 "Poor_Comp_LC")

          for(iter in 1:length(itervec)){
            if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
            True <- readRDS(file.path(x, iter, "True.rds"))

              if(param=="SPR"){
                pred <- calc_SPR(SB0=Report$SB0, R0=exp(Report$beta), SB_t=Report$SB_t_hat, R_t=Report$R_t_hat)
                true <- calc_SPR(SB0=True$SB0, R0=lh_dat$R0, SB_t=True$SB_t, R_t=True$R_t)                
                re[iter,] <- (pred - true)/true
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=10, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,] <- ((Report$F_t/pred) - (True$F_t/true))/(True$F_t/true)
              } 
              if(param=="Depl"){
                pred <- Report$Depl
                true <- True$D_t
                re[iter,] <- (pred - true)/true
              }
            
          }
       } ## end explore_comp


        Outs <- NULL
        Outs$re <- re
        return(Outs)
	}

	ymax <- 3
  	ymin <- -1
	if(length(data_vec)>1){
		par(mfrow=c(nrow, ncol), mar=c(0,0,0,0), omi=c(2,2,1,1))

		if(length(F_scenarios)==1){
			for(r in 1:length(R_scenarios)){
				for(d in 1:length(data_vec)){
          if(simdir=="base"){
   					out <- calc_re(F_choose=F_scenarios, dat_choose=data_vec[d],
   						R_choose=R_scenarios[r], param=param, 
   						simdir=simdir)
          }
          if(simdir=="explore_comp"){
            out <- calc_re(F_choose=F_scenarios, dat_choose=data_vec[d], 
              R_choose=R_scenarios[r], param=param, simdir=simdir)
          }
					if(param=="SPR"){
						color <- "goldenrod"
					}
					if(param=="FFref"){
						color <- "deepskyblue"
					}
          if(param=="Depl"){
            color <- "gray"
          }
					plot_df <- as.data.frame(out$re)
					beanplot(x=c(1:ncol(plot_df)), plot_df, what=c(0,1,1,0),
            			cex=0.5, beanlinewd=3, maxstripline=0.5, maxwidth=0.75, 
            			ylim=c(ymin,ymax), xlim=c(1.5,ncol(plot_df)+1.5), log="",overalline="median", beanlines="median",
            			col=color, xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i",
            			na.rm=TRUE, border=NA, cutmin=min(plot_df, na.rm=TRUE), cutmax=max(plot_df, na.rm=TRUE))
					abline(h=0, lty=2)
					if(r==1) axis(2, at=seq(ymin+1,ymax,by=1), las=2, cex.axis=1.5)
					if(length(R_scenarios)>1) if(r==1 & d==1) mtext(side=3, paste0("R ", R_scenarios[1]), line=1, cex=1.5)
					if(r==2 & d==1) mtext(side=3, paste0("F ", F_scenarios[2]), line=1, cex=1.5)
        			if(r==3 & d==1) mtext(side=3, paste0("F ", F_scenarios[3]), line=1, cex=1.5)

					    if(r==1) mtext(side=2, data_scenario_names[d], line=4, cex=1.5)
        			if(d==length(data_vec)) axis(1, at=seq(5,20,by=5), cex.axis=1.5)
                        mtext("Year", side=1, line=3,  outer=TRUE, cex=1.5)

				}
			}
		}
	}

}