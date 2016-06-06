interval_coverage <- function(Fscen, Rscen, data_scenarios, SigR, 
	itervec, simdir, lh_num, lh_dat, Ncomp_scenarios=FALSE, param,
  sens_param=FALSE, sens_val=FALSE, RecType, OMselex, EMselex){


    if(RecType==0) RecTypeDesc <- "shared"
    if(RecType==1) RecTypeDesc <- "BH"
    if(RecType==2) RecTypeDesc <- "constant"


	if(simdir=="base"){
			
			mat <- rep(0, length(data_scenarios))

	    for(dat in 1:length(data_scenarios)){
    		x <- file.path(test_dir, "base", paste0("LH_", lh_num),
    			paste0("OM", OMselex, "_EM", EMselex), 
    			paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR),
    			paste0("Fdyn_", Fscen), paste0("Rdyn_", Rscen),
    			data_scenarios[dat])
    		for(iter in 1:length(itervec)){
       		  if(file.exists(file.path(x, iter, "True.rds"))) True <- readRDS(file.path(x, iter, "True.rds"))
       		  if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
       		  if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "True.rds"))==FALSE) next
       		  if(file.exists(file.path(x, iter, "Sdreport.rds"))==FALSE) next
       		  if(all(is.na(Sdreport))) next

          if(param=="Depl"){
       		  t <- True$D_t[length(True$D_t)]
       		  s <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
            high <- exp(s[nrow(s),1]+0.67449*s[nrow(s),2])
            low <- exp(s[nrow(s),1]-0.67449*s[nrow(s),2])
          }
          if(param=="FFref"){
              traw <- True$F_t[length(True$F_t)]
              sraw <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]
              rref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                 M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
              tref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                 M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root

              t <- (traw/tref)
              s <- sraw

            high <- exp(s[nrow(s),1]+0.67449*s[nrow(s),2])/rref
            low <- exp(s[nrow(s),1]-0.67449*s[nrow(s),2])/rref

          }

       		  if(is.na(high)|is.na(low)) next
       		  if(t>low & t<high) mat[dat] <- mat[dat]+1
    		}
    		mat[dat] <- mat[dat]/length(itervec)
    	}
    }

    if(simdir=="explore_comp"){

    		mat <- rep(0, length(Ncomp_scenarios))

	    for(dat in 1:length(Ncomp_scenarios)){
    		x <- file.path(test_dir, "explore_comp", paste0("Nyrs_comp_", Ncomp_scenarios[dat]),
    		    paste0("LH_", lh_num),
    			paste0("OM", OMselex, "_EM", EMselex), 
    			paste0("RecType_", RecTypeDesc), paste0("SigR_", SigR),
    			paste0("Fdyn_", Fscen), paste0("Rdyn_", Rscen),
    			"Poor_Comp_LC")
    		for(iter in 1:length(itervec)){
       		  if(file.exists(file.path(x, iter, "True.rds"))) True <- readRDS(file.path(x, iter, "True.rds"))
       		  if(file.exists(file.path(x, iter, "Sdreport.rds"))) Sdreport <- readRDS(file.path(x, iter, "Sdreport.rds"))
       		  if(file.exists(file.path(x, iter, "Report.rds"))) Report <- readRDS(file.path(x, iter, "Report.rds"))
            if(file.exists(file.path(x, iter, "True.rds"))==FALSE) next
       		  if(file.exists(file.path(x, iter, "Sdreport.rds"))==FALSE) next
       		  if(all(is.na(Sdreport))) next


          if(param=="Depl"){
       		  t <- True$D_t[length(True$D_t)]
       		  s <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]
            high <- exp(s[nrow(s),1]+0.67449*s[nrow(s),2])
            low <- exp(s[nrow(s),1]-0.67449*s[nrow(s),2])

          }
          if(param=="FFref"){
              traw <- True$F_t[length(True$F_t)]
              sraw <- summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]
              rref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=Report$Mat_a, W_a=Report$W_a, 
                 M=Report$M, S_a=Report$S_a, R0=exp(Report$beta), ref=0.3)$root
              tref <- uniroot(calc_Fref, lower=0, upper=50, Mat_a=True$Mat_a, W_a=True$W_a, 
                 M=True$M, S_a=True$S_a, R0=True$R0, ref=0.3)$root

              t <- (traw/tref)
              s <- sraw

            high <- exp(s[nrow(s),1]+0.67449*s[nrow(s),2])/rref
            low <- exp(s[nrow(s),1]-0.67449*s[nrow(s),2])/rref

          }
       		  if(is.na(high)|is.na(low)) next
       		  if(t>low & t<high) mat[dat] <- mat[dat]+1
    		}
    		mat[dat] <- mat[dat]/length(itervec)
    	}
    }

    return(mat)

}