EMprocess <- function(dat_scenario, seed, itervec, SigmaR, Fdynamics, Rdynamics, 
	Nyears_comp=FALSE, simdir, lh_dat, lh_num, sens_param=FALSE, sens_val=FALSE, 
  RecType, OMselex, EMselex, est_sigma, rerun=FALSE, runLBSPR=FALSE, REML){
  
	if(Nyears_comp == FALSE) Ncomp <- NULL
	if(Nyears_comp != FALSE) Ncomp <- Nyears_comp

  dat_input <- create_inputs(simdir=simdir, param=sens_param, val=sens_val,
  	lh_dat=lh_dat, RecType=RecType)
  if(sens_param=="SigmaR"){
  	sens_label <- ifelse(dat_input[["SigmaR"]] < SigmaR, "low", 
  		ifelse(dat_input[["SigmaR"]] > SigmaR, "high", FALSE))
  }
  if(sens_param %in% names(lh_dat)){
  	sens_label <- ifelse(dat_input[[sens_param]] < lh_dat[[sens_param]], "low",
  		ifelse(dat_input[[sens_param]] > lh_dat[[sens_param]], "high", FALSE))
  }

  if(RecType==0) RecTypeDesc <- "shared"
  if(RecType==1) RecTypeDesc <- "BH"
  if(RecType==2) RecTypeDesc <- "constant"
  modpath <- setPaths(lh=lh_num, Fdyn=Fdynamics, Rdyn=Rdynamics, iter=NULL, SimDir=simdir,
    SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, iterpath=FALSE, 
    Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc,
    OMselex=OMselex, EMselex=EMselex)
  if(rerun==TRUE){
    unlink(modpath, TRUE)
    modpath <- setPaths(lh=lh_num, Fdyn=Fdynamics, Rdyn=Rdynamics, iter=NULL, SimDir=simdir,
      SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, iterpath=FALSE, 
      Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc,
      OMselex=OMselex, EMselex=EMselex)
  }

    ## settings for each data scenario
    scenario_setup <- NULL
      if(grepl("Rich", dat_scenario)){
      	scenario_setup$Nyears <- 20
      	scenario_setup$comp_sample <- 500
      	scenario_setup$obs_per_yr <- rep(500, scenario_setup$Nyears)
      	scenario_setup$Nyears_comp <- FALSE
      	scenario_setup$alt_yrs <- FALSE
      	scenario_setup$sample <- FALSE
      }

      if(grepl("Moderate_Full", dat_scenario)){
      	scenario_setup$Nyears <- 20
      	scenario_setup$comp_sample <- 500
      	scenario_setup$obs_per_yr <- rep(50, scenario_setup$Nyears)
      	scenario_setup$Nyears_comp <- FALSE
      	scenario_setup$alt_yrs <- FALSE
      	scenario_setup$sample <- FALSE
      }

      if(grepl("Moderate_Sample", dat_scenario)){
      	scenario_setup$Nyears <- 20
      	scenario_setup$comp_sample <- 500
      	scenario_setup$obs_per_yr <- rep(50, scenario_setup$Nyears)
      	scenario_setup$Nyears_comp <- FALSE
      	scenario_setup$alt_yrs <- TRUE
      	scenario_setup$sample <- 0.2
      }

      if(grepl("Poor_Index", dat_scenario)){
      	scenario_setup$Nyears <- 20
      	scenario_setup$comp_sample <- 100
      	scenario_setup$obs_per_yr <- rep(10, scenario_setup$Nyears)
      	scenario_setup$Nyears_comp <- 1
      	scenario_setup$alt_yrs <- FALSE
      	scenario_setup$sample <- FALSE
      }

      if(grepl("Poor_Comp", dat_scenario) | grepl("Poor_Rel", dat_scenario)){
      	scenario_setup$Nyears <- 20
      	scenario_setup$comp_sample <- 100
      	scenario_setup$obs_per_yr <- rep(10, scenario_setup$Nyears)
      	if(Nyears_comp==FALSE) scenario_setup$Nyears_comp <- 1
      	if(Nyears_comp!=FALSE) scenario_setup$Nyears_comp <- Nyears_comp
      	scenario_setup$alt_yrs <- FALSE
      	scenario_setup$sample <- FALSE
      }

  depl_re <- sigr_re <- rep(NA, length(itervec))
 
  set.seed(seed)
  ## generate data
  for(iter in itervec){ ## start write data file

      iterpath <- setPaths(lh=lh_num, Fdyn=Fdynamics, Rdyn=Rdynamics, iter=iter, SimDir=simdir,
        SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, iterpath=TRUE, 
        Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc,
        OMselex=OMselex, EMselex=EMselex)

      DataList <- SimData_LB(Nyears=scenario_setup$Nyears, AgeMax=lh_dat$AgeMax,
        M=lh_dat$M, F1=lh_dat$F1, h=lh_dat$h, S_a=lh_dat$S_a, qcoef=lh_dat$qcoef, 
        Frate=lh_dat$Frate, Fequil=lh_dat$Fequil, 
        SigmaF=lh_dat$SigmaF,SigmaR=SigmaR, Fdynamics=Fdynamics, L_a=lh_dat$L_a, W_a=lh_dat$W_a,
        Rdynamics=Rdynamics, R0=lh_dat$R0, Fmax=lh_dat$Fmax, CVlen=lh_dat$CVlen, 
        Mat_a=lh_dat$Mat_a, Amat=lh_dat$Amat, linf=lh_dat$linf, vbk=lh_dat$vbk,
        mids=lh_dat$mids, highs=lh_dat$highs, lows=lh_dat$lows, Nyears_comp=scenario_setup$Nyears_comp,
        comp_sample=scenario_setup$comp_sample, alt_yrs=scenario_setup$alt_yrs, 
        sample=scenario_setup$sample, nburn=20) 

       saveRDS(DataList, file.path(iterpath, "True.rds"))
       rm(DataList)
       rm(iterpath)
  } ## end write data file

  dir.create(file.path(modpath, "figs"), showWarnings=FALSE)


  for(iter in itervec){

  	iterpath <- file.path(modpath, iter)

    DataList <- readRDS(file.path(iterpath, "True.rds"))

    if(runLBSPR==TRUE & grepl("Poor_Comp", dat_scenario)){
        if(is.null(Ncomp)){
          lfreq <- data.frame("bin"=1:length(DataList$LF), "obs"=DataList$LF[length(DataList$LF)])
          lendat <- unlist(sapply(1:nrow(lfreq), function(x) rep(lfreq[x,"bin"], lfreq[x,"obs"])))
          lbspr_res <- LBSPR_assessment(lh_list=lh_dat, lendat=lendat, setSPR=0.3)
          saveRDS(lbspr_res, file.path(iterpath, "LBSPR_results.rds"))
        }
        if(is.null(Ncomp)==FALSE){
          for(nn in 1:Ncomp){
            if(Ncomp>1) lfreq <- data.frame("bin"=1:ncol(DataList$LF), "obs"=DataList$LF[(nrow(DataList$LF)-(Ncomp-nn)),])
            if(Ncomp==1) lfreq <- data.frame("bin"=1:length(DataList$LF), "obs"=DataList$LF[length(DataList$LF)])
            lendat <- unlist(sapply(1:nrow(lfreq), function(x) rep(lfreq[x,"bin"], lfreq[x,"obs"])))
            lbspr_res <- LBSPR_assessment(lh_list=lh_dat, lendat=lendat, setSPR=0.3)
            saveRDS(lbspr_res, file.path(iterpath, paste0("LBSPR_results_", nn, ".rds")))
          }
        }

    }

    if(file.exists(file.path(iterpath, "Report.rds"))) next


    if(RecType!=2) RecDev_biasadj1 <- rep(1, scenario_setup$Nyears)
    if(RecType==2) RecDev_biasadj1 <- rep(0, scenario_setup$Nyears)
      TmbList1 <- FormatInput_LB(Nyears=scenario_setup$Nyears, Nlenbins=length(dat_input$mids), 
        DataList=DataList, linf=dat_input$linf, vbk=dat_input$vbk, t0=dat_input$t0, M=dat_input$M, AgeMax=dat_input$AgeMax,
        lbhighs=dat_input$highs, lbmids=dat_input$mids, Mat_a=dat_input$Mat_a, 
        lwa=dat_input$lwa, lwb=dat_input$lwb, log_sigma_C=dat_input$log_sigma_C, 
        log_sigma_I=dat_input$log_sigma_I, log_CV_L=dat_input$log_CV_L, F1=dat_input$F1,
        SigmaR=SigmaR, qcoef=dat_input$qcoef, R0=dat_input$R0, S50=dat_input$S50, S95=dat_input$S95, 
        version=Version, model=dat_scenario, RecDev_biasadj=RecDev_biasadj1, SigmaF=dat_input$SigmaF,
        Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), obs_per_yr=scenario_setup$obs_per_yr, RecType=RecType, h=dat_input$h,
        SelexTypeDesc=EMselex, est_sigma=est_sigma, REML=REML)   

      saveRDS(TmbList1, file.path(iterpath, "Inputs1.rds"))
      # TmbList1 <- readRDS(file.path(iterpath, "Inputs1.rds"))


      dyn.load(paste0(run_exe, "\\", dynlib(Version)))
      obj1 <- MakeADFun(data=TmbList1[["Data"]], parameters=TmbList1[["Parameters"]],
                        random=TmbList1[["Random"]], map=TmbList1[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE)       

      ## Settings
      obj1$env$inner.control <- c(obj1$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
      obj1$hessian <- FALSE 
        Upr = rep(Inf, length(obj1$par))
          Upr[match("log_sigma_R",names(obj1$par))] = log(2)
          Upr[match("S50", names(obj1$par))] = dat_input$AgeMax
          Upr[which(names(obj1$par)=="log_F_t_input")] = log(2)
          Upr[match("log_F_sd", names(obj1$par))] <- log(2)
          Lwr = rep(-Inf, length(obj1$par))
          Lwr[which(names(obj1$par)=="log_F_t_input")] = log(0.001) 
          Lwr[which(names(obj1$par)=="log_sigma_R")] = log(0.001)
          Lwr[which(names(obj1$par)=="log_F_sd")] <- log(0.001)   

      ## Run optimizer
      opt1 <- tryCatch( nlminb( start=obj1$par, objective=obj1$fn, gradient=obj1$gr, upper=Upr, lower=Lwr,
        control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)
      if(all(is.na(opt1))==FALSE){
        opt1[["final_gradient"]] =  obj1$gr( opt1$par )      
        df1 <- data.frame(opt1$final_gradient, names(obj1$par), opt1$par)  
      }    

      ## Standard errors
      Report1 = obj1$report()  
      Sdreport1 = tryCatch( sdreport(obj1), error=function(x) NA )
      ParList <- obj1$env$parList( obj1$env$last.par.best )
      # dyn.unload( paste0(run_exe,"\\", dynlib(Version)) )      

    if(all(is.na(Sdreport1))) RecDev_biasadj <- rep(1, scenario_setup$Nyears)
	  if(all(is.na(Sdreport1))==FALSE){
        SD <- summary(Sdreport1)
        if(RecType!=2) RecDev_biasadj <- 1 - SD[which(rownames(SD)=="Nu_input"), "Std. Error"]^2 / Report1$sigma_R^2    
        if(RecType==2) RecDev_biasadj <- rep(0, scenario_setup$Nyears)
	  }

      TmbList <- FormatInput_LB(Nyears=scenario_setup$Nyears, Nlenbins=length(dat_input$mids), 
        DataList=DataList, linf=dat_input$linf, vbk=dat_input$vbk, t0=dat_input$t0, M=dat_input$M, AgeMax=dat_input$AgeMax,
        lbhighs=dat_input$highs, lbmids=dat_input$mids, Mat_a=dat_input$Mat_a, 
        lwa=dat_input$lwa, lwb=dat_input$lwb, log_sigma_C=dat_input$log_sigma_C, log_sigma_I=dat_input$log_sigma_I, log_CV_L=dat_input$log_CV_L, F1=dat_input$F1,
        SigmaR=SigmaR, qcoef=dat_input$qcoef, R0=dat_input$R0, S50=dat_input$S50, S95=dat_input$S95, 
        version=Version, model=dat_scenario, RecDev_biasadj=RecDev_biasadj, SigmaF=dat_input$SigmaF,
        Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), obs_per_yr=scenario_setup$obs_per_yr, RecType=RecType, h=dat_input$h,
        SelexTypeDesc=EMselex, est_sigma=est_sigma, REML=REML)

      saveRDS(TmbList, file.path(iterpath, "Inputs2.rds"))

      dyn.load(paste0(run_exe, "\\", dynlib(Version)))
      obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE)       

      InitVal <- obj$fn(obj$par)        

      ## Settings
      obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
      obj$hessian <- FALSE 
        Upr = rep(Inf, length(obj$par))
          Upr[match("log_sigma_R",names(obj$par))] = log(2)
          Upr[match("S50", names(obj$par))] = dat_input$AgeMax
          Upr[which(names(obj$par)=="log_F_t_input")] = log(2)
          Upr[match("log_F_sd", names(obj$par))] <- log(2)
          Lwr = rep(-Inf, length(obj$par))
          Lwr[which(names(obj$par)=="log_F_t_input")] = log(0.001) 
          Lwr[which(names(obj$par)=="log_sigma_R")] = log(0.001)
          Lwr[which(names(obj$par)=="log_F_sd")] <- log(0.001)
    

      ## Run optimizer
      opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr,
        control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)
      for(i in 1:5){
      	if(all(is.na(opt))){
      		obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
      	} else{
      		break
      	}
      }
      if(all(is.na(opt))){
      	write("NAs final gradient", file.path(iterpath, "NAs_final_gradient.txt"))
		    next
	    }
      if(all(is.na(opt))==FALSE){

      	opt[["final_gradient"]] = obj$gr( opt$par )       
      	df <- data.frame(opt$final_gradient, names(obj$par), opt$par)      

      	for(i in 1:5){
          if(abs(min(opt[["final_gradient"]]))>0.01){
            obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
            if(all(is.na(opt))==FALSE) opt[["final_gradient"]] <- obj$gr(opt$par) 
            if(all(is.na(opt))) break
          } else{
            break
          } 
        }
        for(i in 1:5){
      	  if(all(is.na(opt))){
      		obj <- MakeADFun(data=TmbList[["Data"]], parameters=ParList,
                        random=TmbList[["Random"]], map=TmbList[["Map"]], 
                        inner.control=list(maxit=1e3), hessian=FALSE) 
            opt <- tryCatch( nlminb( start= obj$env$last.par.best[-obj$env$random] + rnorm(length(obj$par),0,0.2), 
              objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, 
              control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10)), error=function(e) NA)
      	  } else{
      		break
      	  }
        }
        if(all(is.na(opt))){
      	  write("NAs final gradient", file.path(iterpath, "NAs_final_gradient.txt"))
		  next
	    } 
        if(abs(min(opt[["final_gradient"]]))>0.01) write(opt[["final_gradient"]], file.path(iterpath, "high_final_gradient.txt"))
      }      

      ## Standard errors
      Report = tryCatch( obj$report(), error=function(x) NA)
  	  saveRDS(Report, file.path(iterpath, "Report.rds"))

      Sdreport = tryCatch( sdreport(obj), error=function(x) NA )
      saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))

      	depl_re[iter-itervec[1]+1] <- (Report$SB_t_hat[length(Report$SB_t_hat)]/Report$SB0 - DataList$D_t[length(DataList$D_t)])/DataList$D_t[length(DataList$D_t)]
        sigr_re[iter-itervec[1]+1] <- (exp(Report$log_sigma_R) - DataList$SigmaR)/DataList$SigmaR

      if(iter==1 | iter %% 10 == 0 & all(is.na(Sdreport))==FALSE){
      	figname <- paste0(iter,"_Model_fit.png")
        png(file=file.path(modpath, "figs", figname), width=9, height=6, res=200, units="in")

      	par(mfrow=c(3,3), mar=c(3,3,2,0))
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1*InputMat[,2]), rev(exp(InputMat[,1]+1*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1*InputMat[,2], rev(InputMat[,1]+1*InputMat[,2])))
        } 
        ## Abundance
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$SB_t, "Est"=Report$SB_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Biomass", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Index
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$N_t, "Est"=Report$N_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Abundance", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Average Length
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$L_t, "Est"=Report$L_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Average Length", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Recruitment
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$R_t, "Est"=Report$R_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)    
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Recruitment", lwd=2)
        if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Fishing Mortality
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$F_t, "Est"=Report$F_t)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Fishing Mortality", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Depletion
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=DataList$D_t, "Est"=Report$Depl)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Depletion", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Catch
        tcatch <- DataList$C_t
        pcatch <- rep(NA, length(DataList$SB_t))
          names(pcatch) <- 1:length(DataList$SB_t)
        pcatch[which(names(pcatch) %in% names(DataList$C_t))] <- tcatch
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=pcatch, "Est"=Report$C_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Catch", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
        ## Index
        tindex <- DataList$I_t
        pindex <- rep(NA, length(DataList$SB_t))
          names(pindex) <- 1:length(DataList$SB_t)
        pindex[which(names(pindex) %in% names(DataList$I_t))] <- tindex
        Mat <- cbind("Year"=1:length(DataList$SB_t), "True"=pindex, "Est"=Report$I_t_hat)
        ymax <- ifelse(max(Mat[,c("True", "Est")], na.rm=TRUE)==Inf, 1, max(Mat[,c("True", "Est")], na.rm=TRUE)*1.5)
        matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), main="Index", lwd=2)
        if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)

        # text(x=15, y=20000, paste0("sigmaR RE = ", round((Report$sigma_R - SigmaR)/SigmaR, 3)))
        # text(x=15, y=10000, paste0("depl RE = ", round((Report$SB_t_hat[length(Report$SB_t_hat)]/Report$SB0 - DataList$D_t[length(DataList$D_t)])/DataList$D_t[length(DataList$D_t)], 3)))

        dev.off()

        dyn.unload( paste0(run_exe,"\\", dynlib(Version)) ) 

        write.csv(df, file.path(modpath, "df.csv"))

      } 

  		rm(Report)
  		rm(Sdreport)
  		rm(TmbList)
  		rm(TmbList1)
  		rm(opt)
  		rm(opt1)
  		rm(obj)
  		rm(obj1)
  		rm(df)
  		rm(df1)
  		rm(InitVal)
  		rm(ParList)

    }

	Outs <- NULL
	Outs$depl_re <- depl_re
	Outs$sigr_re <- sigr_re
	return(Outs)


}
