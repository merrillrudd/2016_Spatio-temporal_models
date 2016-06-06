REplots <- function(F_scenarios, R_scenarios, data_scenarios, SigR, 
	itervec, param, simdir, lh_num, lh_dat, Ncomp_scenarios=FALSE,
  sens_param=FALSE, sens_val=FALSE, RecType, OMselex, EMselex, retype=1){

	ncol <- length(F_scenarios)
	nrow <- length(R_scenarios)

	par(mfcol=c(ncol, nrow), mar=c(0,0,0,0), omi=c(2,1,1,1))

	calc_re <- function(F_choose, R_choose, param, simdir, RecType){
     if(simdir=="base"){
     	re <- matrix(NA, nrow=length(itervec), ncol=length(data_scenarios))
     	NAsSD <- HighGradient <- matrix(0, nrow=length(itervec), ncol=length(data_scenarios))

      if(RecType==0) RecTypeDesc <- "shared"
      if(RecType==1) RecTypeDesc <- "BH"
      if(RecType==2) RecTypeDesc <- "constant"

        for(dat in 1:length(data_scenarios)){
              x <- file.path(test_dir, "base", paste0("LH_", lh_num), paste0("OM", OMselex, "_EM", EMselex),
                paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR), paste0("Fdyn_", F_choose), 
                paste0("Rdyn_", R_choose), data_scenarios[dat])
          for(iter in 1:length(itervec)){
            if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))==FALSE) NAsSD[iter,dat] <- 1
            if(file.exists(file.path(x, iter, "high_final_gradient.txt"))){
            	HighGradient[iter,dat] <- 1
            	next
            }
            if(NAsSD[iter,dat]==1) next
            True <- readRDS(file.path(x, iter, "True.rds"))

            if(retype==1){
              if(param=="SigmaR"){
              	if(identical(Report$sigma_R, True$SigmaR)) next
              	re[iter,dat] <- (exp(Report$log_sigma_R) - True$SigmaR)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- (Report$Depl[length(Report$Depl)] - True$D_t[length(True$D_t)])/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- ((pred-true)/true)
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- ((Report$F_t[length(Report$F_t)]/pred) - (True$F_t[length(True$F_t)]/true))/(True$F_t[length(True$F_t)]/true)
              } 
            }
            if(retype==2){
              if(param=="SigmaR"){
                if(identical(Report$sigma_R, True$SigmaR)) next
                re[iter,dat] <- exp(Report$log_sigma_R)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- Report$Depl[length(Report$Depl)]/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- pred/true
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- (Report$F_t[length(Report$F_t)]/pred)/(True$F_t[length(True$F_t)]/true)
              } 
            }
          }
        }
      } ## end base

      if(simdir=="explore_comp"){
      	if(length(data_scenarios) > 1) stop("Can only take 1 data scenario")
        
        re <- matrix(NA, nrow=length(itervec), ncol=length(Ncomp_scenarios))
     	NAsSD <- HighGradient <- matrix(0, nrow=length(itervec), ncol=length(Ncomp_scenarios))

      if(RecType==0) RecTypeDesc <- "shared"
      if(RecType==1) RecTypeDesc <- "BH"
      if(RecType==2) RecTypeDesc <- "constant"


        for(dat in 1:length(Ncomp_scenarios)){

          	x <- file.path(test_dir, "explore_comp", paste0("Nyrs_comp_", Ncomp_scenarios[dat]), 
          	     paste0("LH_", lh_num),  paste0("OM", OMselex, "_EM", EMselex), paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR), 
                 paste0("Fdyn_", F_choose), paste0("Rdyn_", R_choose),
          	     data_scenarios)

          for(iter in 1:length(itervec)){
            if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))==FALSE) NAsSD[iter,dat] <- 1
            if(file.exists(file.path(x, iter, "high_final_gradient.txt"))){
            	HighGradient[iter,dat] <- 1
            	next
            }
            if(NAsSD[iter,dat]==1) next
            True <- readRDS(file.path(x, iter, "True.rds"))

            if(retype==1){
              if(param=="SigmaR"){
                if(identical(Report$sigma_R, True$SigmaR)) next
                re[iter,dat] <- (exp(Report$log_sigma_R) - True$SigmaR)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- (Report$Depl[length(Report$Depl)] - True$D_t[length(True$D_t)])/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- ((pred-true)/true)
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- ((Report$F_t[length(Report$F_t)]/pred) - (True$F_t[length(True$F_t)]/true))/(True$F_t[length(True$F_t)]/true)
              } 
            }
            if(retype==2){
              if(param=="SigmaR"){
                if(identical(Report$sigma_R, True$SigmaR)) next
                re[iter,dat] <- exp(Report$log_sigma_R)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- Report$Depl[length(Report$Depl)]/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- pred/true
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- (Report$F_t[length(Report$F_t)]/pred)/(True$F_t[length(True$F_t)]/true)
              } 
            } 
          }
        }
      } ## end explore comp

     if(simdir=="sensitivity" & all(Ncomp_scenarios==FALSE)){
      re <- matrix(NA, nrow=length(itervec), ncol=length(data_scenarios))
      NAsSD <- HighGradient <- matrix(0, nrow=length(itervec), ncol=length(data_scenarios))

      if(RecType==0) RecTypeDesc <- "shared"
      if(RecType==1) RecTypeDesc <- "BH"
      if(RecType==2) RecTypeDesc <- "constant"


        for(dat in 1:length(data_scenarios)){
          x <- file.path(test_dir, "sensitivity", paste0(sens_param, "_", sens_val), 
            paste0("LH_", lh_num), paste0("OM", OMselex, "_EM", EMselex),  
            paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR), 
            paste0("Fdyn_", F_choose), paste0("Rdyn_", R_choose),
            data_scenarios[dat])
          for(iter in 1:length(itervec)){
            if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
            if(file.exists(file.path(x, iter, "Sdreport.rds"))==FALSE) NAsSD[iter,dat] <- 1
            if(file.exists(file.path(x, iter, "high_final_gradient.txt"))){
              HighGradient[iter,dat] <- 1
              next
            }
            if(NAsSD[iter,dat]==1) next
            True <- readRDS(file.path(x, iter, "True.rds"))

            if(retype==1){
              if(param=="SigmaR"){
                if(identical(Report$sigma_R, True$SigmaR)) next
                re[iter,dat] <- (exp(Report$log_sigma_R) - True$SigmaR)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- (Report$Depl[length(Report$Depl)] - True$D_t[length(True$D_t)])/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- ((pred-true)/true)
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- ((Report$F_t[length(Report$F_t)]/pred) - (True$F_t[length(True$F_t)]/true))/(True$F_t[length(True$F_t)]/true)
              } 
            }
            if(retype==2){
              if(param=="SigmaR"){
                if(identical(Report$sigma_R, True$SigmaR)) next
                re[iter,dat] <- exp(Report$log_sigma_R)/True$SigmaR
              }
              if(param=="Depl") re[iter,dat] <- Report$Depl[length(Report$Depl)]/True$D_t[length(True$D_t)]
              if(param=="SPR"){
                pred <- calc_ref(Mat_a=Report$Mat_a, W_a=Report$W_a, M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), F=Report$F_t[length(Report$F_t)], ref=FALSE)
                true <- calc_ref(Mat_a=True$Mat_a, W_a=True$W_a, M=True$M, S_a=True$S_a, R0=exp(True$beta), F=True$F_t[length(True$F_t)], ref=FALSE)
                re[iter,dat] <- pred/true
              }
              if(param=="FFref"){
                pred <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                  M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
                true <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                  M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root
                re[iter,dat] <- (Report$F_t[length(Report$F_t)]/pred)/(True$F_t[length(True$F_t)]/true)
              } 
            } 
          }
        }
      } ## end sensitivity


        Outs <- NULL
        Outs$re <- re
        Outs$NAsSD <- NAsSD
        Outs$HighGradient <- HighGradient
        return(Outs)
	}

if(retype==1){
	ymax <- 3
  ymin <- -1
}
if(retype==2){
  ymin=0
  ymax=4
}
      if(length(param)==1){
        if(param=="Depl") color <- "gray"
        if(param=="SigmaR") color <- "goldenrod"
        if(param=="FFref") color <- "deepskyblue"
        if(param=="SPR") color <- "tomato"
      }
      if(length(param)>1){
        color <- list()
        color[[1]] <- "black"
        if("Depl" %in% param) color[[which(param=="Depl")+1]] <- "#AAAAAA"
        if("SigmaR" %in% param) color[[which(param=="SigmaR")+1]] <- "#AA0000"
        if("FFref" %in% param) color[[which(param=="FFref")+1]] <- "#00AA00"
        if("SPR" %in% param) color[[which(param=="SPR")+1]] <- "#AA00AA"
      }

	for(f in 1:length(F_scenarios)){
		for(r in 1:length(R_scenarios)){

      if(length(param)==1){
        out <- calc_re(F_choose=F_scenarios[f], R_choose=R_scenarios[r], param=param, simdir=simdir, RecType=RecType)
        plot_df <- as.data.frame(out$re)
      }
      if(length(param)>1){
        for(p in 1:length(param)){
          out <- calc_re(F_choose=F_scenarios[f], R_choose=R_scenarios[r], param=param[p], simdir=simdir, RecType=RecType)
          if(p==1) plot_df <- as.data.frame(out$re)
          if(p>1) plot_df <- cbind.data.frame(plot_df, as.data.frame(out$re))
        }
      }

      ## remove NAs
      for(p in 1:length(param)){
        index <- which(is.na(plot_df[,p])==FALSE)
        plot_df <- plot_df[index,]
      }

			tryCatch(beanplot(x=c(1:ncol(plot_df)), plot_df, what=c(0,1,1,0),
            cex=0.5, beanlinewd=3, maxstripline=0.5, maxwidth=0.75, 
            ylim=c(ymin,ymax), xlim=c(1.5,ncol(plot_df)+1.5), log="",overalline="median", beanlines="median",
            col=color, xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i",
            na.rm=TRUE, border=NA, cutmin=min(plot_df, na.rm=TRUE), cutmax=max(plot_df, na.rm=TRUE)), error=function(e) print("NAs"))

			# boxplot(out$re, col=color,
  	# 			ylim=c(-1,ymax), xaxt="n", yaxt="n")
			if(retype==1) abline(h=0, lty=2)
      if(retype==2) abline(h=1, lty=2)
			# for(i in 1:ncol(out$re)){
			# 	text(x=i, y=ymax-0.1, round(sum(out$NAsSD[,i])/length(itervec),2), col="red", cex=1.3)
			# 	text(x=i, y=ymax-0.5, round(sum(out$HighGradient[,i])/length(itervec),2), col="blue", cex=1.3)
			# }
			if(f==1) axis(2, at=seq(ymin+1,ymax,by=1), las=2, cex.axis=1.5)
			if(r==3){
        labels <- rep(NA, length(data_scenarios))
        if(length(data_scenarios)>1){
          if("Rich_LC" %in% data_scenarios) labels[which(data_scenarios=="Rich_LC")] <- "Rich"
          if("Moderate_Full_LC" %in% data_scenarios) labels[which(data_scenarios=="Moderate_Full_LC")] <- "Mod\nFull"
          if("Moderate_Sample_LC" %in% data_scenarios) labels[which(data_scenarios=="Moderate_Sample_LC")] <- "Mod\nSample"
          if("Poor_Index_LC" %in% data_scenarios) labels[which(data_scenarios=="Poor_Index_LC")] <- "Poor\nIndex"
          if("Poor_Comp_LC" %in% data_scenarios) labels[which(data_scenarios=="Poor_Comp_LC")] <- "Poor\nComp"
          if("Poor_Rel_LC" %in% data_scenarios) labels[which(data_scenarios=="Poor_Rel_LC")] <- "Poor\nRel"
        }
        if(length(param)>1){
          labels <- param
        }

				if(simdir %in% c("base", "sensitivity")) text(x=2:(length(labels)+1), y=ymin-1, xpd=NA, cex=1.5, labels)
				if(simdir=="explore_comp") axis(1, at=c(1:length(Ncomp_scenarios))+1, cex.axis=1.5, labels=Ncomp_scenarios)
			}
			if(f==1 & r==1) mtext(side=3, paste0("F ", F_scenarios[1]), line=1, cex=1.5)
			if(f==1 & r==1) mtext(side=2, paste0("R ", R_scenarios[1]), line=4, cex=1.5)
			if(f==1 & r==2) mtext(side=2, paste0("R ", R_scenarios[2]), line=4, cex=1.5)
			if(f==1 & r==3) mtext(side=2, paste0("R ", R_scenarios[3]), line=4, cex=1.5)
			if(f==2 & r==1) mtext(side=3, paste0("F ", F_scenarios[2]), line=1, cex=1.5)
			if(f==3 & r==1) mtext(side=3, paste0("F ", F_scenarios[3]), line=1, cex=1.5)
		}
	}
	if(simdir %in% c("base", "sensitivity")) mtext(side=1, outer=TRUE, "Data Scenario", line=8)
  if(simdir=="explore_comp") mtext(side=1, outer=TRUE, "Number of Years of Composition Data", line=3)
	mtext(side=2, outer=TRUE, "Relative Error", line=2)


}