runModel <- function(modpath, itervec, data_input, Fdynamics, RecType=0, FType=0, LType=1, site=1){

Fdyn_dir <- file.path(modpath, paste0("F_", Fdynamics))

## run model 
for(iter in itervec){

	iterpath <- file.path(Fdyn_dir, iter)
	DataList <- readRDS(file.path(iterpath, "True.rds"))

   RecDev_biasadj <- rep(0, Nyears)
   TmbList <- FormatInput_LB(Nyears=Nyears, DataList=DataList, linf=data_input$linf, vbk=data_input$vbk, t0=data_input$t0, M=data_input$M, AgeMax=data_input$AgeMax, lbhighs=data_input$highs, lbmids=data_input$mids, Mat_a=data_input$Mat_a, lwa=data_input$lwa, lwb=data_input$lwb, log_sigma_C=data_input$log_sigma_C, log_sigma_I=data_input$log_sigma_I, log_CV_L=data_input$log_CV_L, F1=data_input$F1, SigmaR=SigmaR, qcoef=data_input$qcoef, R0=data_input$R0, S50=data_input$S50, S95=data_input$S95, version=Version, model=dat_scenario, RecDev_biasadj=RecDev_biasadj,SigmaF=data_input$SigmaF, Fpen=1, Dpen=0, Dprior=c(0.8, 0.2), obs_per_yr=rep(n_i, Nyears), RecType=RecType, FType=FType, LType=LType, h=data_input$h, SelexTypeDesc="asymptotic", est_sigma=est_sigma, REML=FALSE, site=site)
   saveRDS(TmbList, file.path(iterpath, "Inputs.rds"))   

   dyn.load(paste0(run_exe, "\\", dynlib(Version)))   

   obj <- MakeADFun(data=TmbList[["Data"]], parameters=TmbList[["Parameters"]], random=TmbList[["Random"]], map=TmbList[["Map"]],inner.control=list(maxit=1e3), hessian=FALSE)          

   ## Settings
   obj$env$inner.control <- c(obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
   obj$hessian <- FALSE 
     Upr = rep(Inf, length(obj$par))
     Upr[match("log_sigma_R",names(obj$par))] = log(2)
     Upr[match("S50", names(obj$par))] = data_input$AgeMax
     Upr[which(names(obj$par)=="log_F_t_input")] = log(2)
     Upr[match("log_F_sd", names(obj$par))] <- log(2)
     Lwr = rep(-Inf, length(obj$par))
     Lwr[which(names(obj$par)=="log_F_t_input")] = log(0.001) 
     Lwr[which(names(obj$par)=="log_sigma_R")] = log(0.001)
     Lwr[which(names(obj$par)=="log_F_sd")] <- log(0.001)      

   ## Run optimizer
   opt <- tryCatch( nlminb( start=obj$par, objective=obj$fn, gradient=obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) ), error=function(e) NA)   

   if(all(is.na(opt))==FALSE){
       opt[["final_gradient"]] =  obj$gr( opt$par )      
       df <- data.frame(opt$final_gradient, names(obj$par), opt$par)  
   }
   
   ParList <- obj$env$parList( obj$env$last.par.best )
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

        df <- data.frame(opt$final_gradient, names(obj$par), opt$par)      
            

   Report <- obj$report()
   Sdreport <- sdreport(obj)   

   saveRDS(Report, file.path(iterpath, "Report.rds"))
   saveRDS(Sdreport, file.path(iterpath, "Sdreport.rds"))

   rm(Report)
   rm(Sdreport)
   rm(obj)
   rm(opt)
   rm(DataList)
   rm(iterpath)
   rm(TmbList)

}

return(paste0(max(itervec), " iterates run in ", modpath))


}