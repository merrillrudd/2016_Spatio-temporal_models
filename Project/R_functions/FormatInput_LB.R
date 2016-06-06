FormatInput_LB <- function(Nyears, DataList, linf, vbk, t0, M, AgeMax,
	lbhighs, lbmids, Mat_a, lwa, lwb, log_sigma_C, log_sigma_I, log_CV_L, F1, SigmaR, 
	qcoef, R0, S50, S95, version, model, RecDev_biasadj, site,
    Fpen, Dpen, Dprior, obs_per_yr, SigmaF, RecType, FType, LType, h, SelexTypeDesc, est_sigma, REML){

        if(grepl("lb_statespace_v6", version)){
        if(grepl("Poor_Comp", model) & grepl("LC", model)){
            if(is.matrix(DataList$LF) | is.array(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=0, 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                obs_per_yr=obs_per_yr, FType=FType,
                site=site, n_s=length(site),
                I_t=as.vector(0), C_t=as.vector(0), 
                ML_t=as.array(0), LF=LF,
                rel_c=0, rel_i=0,
                vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Comp", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=0, 
                n_lc=0,
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                obs_per_yr=obs_per_yr, FType=FType,
                I_t=as.vector(0), C_t=as.vector(0), 
                ML_t=DataList$ML_t, LF=as.array(0),
                rel_c=0, rel_i=0,
                vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        Parameters <- list(log_F_sd=log(SigmaF), log_F_t_input=log(rep(F1,Nyears)),log_q_I=log(qcoef), beta=log(R0), 
            log_sigma_R=log(SigmaR), S50=S50, S95=S95, log_sigma_C=log_sigma_C, log_sigma_I=log_sigma_I, 
            log_CV_L=log_CV_L,log_sigma_linf=log_CV_L, log_linf=log(linf), Nu_input=rep(0,Nyears),Eps_input=rep(0,length(site)))

        Map = list()
            Map[["log_F_sd"]] <- NA
            Map[["log_F_sd"]] <- factor(Map[["log_F_sd"]])
            Map[["S95"]] <- NA
            Map[["S95"]] <- factor(Map[["S95"]])
            Map[["log_linf"]] <- NA
            Map[["log_linf"]] <- factor(Map[["log_linf"]])
            Map[["log_sigma_linf"]] <- NA
            Map[["log_sigma_linf"]] <- factor(Map[["log_sigma_linf"]])



        if(grepl("Poor_Comp", model)){
            Map[["log_q_I"]] <- NA
            Map[["log_q_I"]] <- factor(Map[["log_q_I"]])

            if(n_lc==1){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][1:10] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
            }
            if(n_lc>1 & n_lc<=10){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][1:10] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
            }

            Map[["log_sigma_R"]] <- NA
            Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

            Map[["beta"]] <- NA
            Map[["beta"]] <- factor(Map[["beta"]])
        } 
        ## turn off annual F estimates - only last value
        #### SHOULD CHANGE THIS TO MEAN-LENGTH MORTALITY ESTIMATOR
        if(FType==1){
            Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
            Map[["log_F_t_input"]] <- NA
            Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
        }
        if(RecType==0 & LType==0) Random <- c("Nu_input", "Eps_input")
        if(RecType==0 & LType==1){
            Random <- c("Nu_input")
            Map[["Eps_input"]] <- 1:length(Parameters[["Eps_input"]])
            Map[["Eps_input"]] <- factor(Map[["Eps_input"]])
        }
        if(RecType==1 & LType==0){
            Random <- c("Eps_input")
            Map[["Nu_input"]] <- 1:length(Parameters[["Nu_input"]])
            Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
        }
        if(RecType==1 & LType==1){
            Random <- NULL
            Map[["Nu_input"]] <- 1:length(Parameters[["Nu_input"]])
            Map[["Nu_input"]] <- factor(Map[["Nu_input"]])
            Map[["Eps_input"]] <- 1:length(Parameters[["Eps_input"]])
            Map[["Eps_input"]] <- factor(Map[["Eps_input"]])

        }


    }
  

        if(grepl("lb_statespace_v5", version)){
        if(grepl("Poor", model)==FALSE & grepl("LC", model)){
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=length(DataList$I_t), 
                n_lc=nrow(DataList$LF),
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.numeric(rownames(DataList$LF)),
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=DataList$C_t, 
                ML_t=as.vector(0), LF=DataList$LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor", model)==FALSE & grepl("ML", model)){
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=length(DataList$C_t),
                n_i=length(DataList$I_t), 
                n_lc=0,
                n_ml=nrow(DataList$LF),
                T_yrs=1:Nyears, C_yrs=as.numeric(names(DataList$C_t)),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=as.numeric(rownames(DataList$LF)),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=DataList$C_t, 
                ML_t=rowMeans(DataList$LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Index", model) & grepl("LC", model)){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=length(DataList$I_t), 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0, 
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=as.vector(0), 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Index", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=length(DataList$I_t), 
                n_lc=0,
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.numeric(names(DataList$I_t)),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=DataList$I_t, C_t=as.vector(0), 
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Comp", model) & grepl("LC", model)){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=0, 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=as.vector(0), 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen,  Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Comp", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=0,
                n_i=0, 
                n_lc=0,
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=as.vector(0),
                I_yrs=as.vector(0),
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                rel_c=0, rel_i=0,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=as.vector(0), C_t=as.vector(0), 
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Rel", model) & grepl("LC", model)){
            if(is.matrix(DataList$LF)){
                n_lc <- nrow(DataList$LF)
                LC_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_lc <- 1
                LC_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            C_yrs <- c(Nyears-19)
            I_yrs <- c(Nyears-19)
            I_t <- c(DataList$I_t[Nyears-19])/DataList$I_t[Nyears]
            C_t <- c(DataList$C_t[Nyears-19])/DataList$C_t[Nyears]

            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=1,
                n_i=1, 
                n_lc=n_lc,
                n_ml=0,
                T_yrs=1:Nyears, C_yrs=C_yrs,
                I_yrs=I_yrs,
                LC_yrs=LC_yrs,
                ML_yrs=as.vector(0),
                rel_c=1, rel_i=1,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=I_t, C_t=C_t, 
                ML_t=as.vector(0), LF=LF,
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }

        if(grepl("Poor_Rel", model) & grepl("ML", model)){
            if(is.matrix(DataList$LF)){
                n_ml <- nrow(DataList$LF)
                ML_yrs <- as.numeric(rownames(DataList$LF))
                LF <- as.matrix(DataList$LF)
            }
            if(is.vector(DataList$LF)){
                n_ml <- 1
                ML_yrs <- Nyears
                LF <- t(as.matrix(DataList$LF))
            }
            C_yrs <- c(Nyears-19)
            I_yrs <- c(Nyears-19)
            I_t <- c(DataList$I_t[Nyears-19])/DataList$I_t[Nyears]
            C_t <- c(DataList$C_t[Nyears-19])/DataList$C_t[Nyears]

            Data <- list(n_t=Nyears, n_lb=ncol(DataList$LF), 
                n_c=1,
                n_i=1, 
                n_lc=0,
                n_ml=n_ml,
                T_yrs=1:Nyears, C_yrs=C_yrs,
                I_yrs=I_yrs,
                LC_yrs=as.vector(0),
                ML_yrs=ML_yrs,
                rel_c=1, rel_i=1,
                obs_per_yr=obs_per_yr, RecType=RecType,
                I_t=I_t, C_t=C_t, 
                ML_t=rowMeans(LF), LF=as.matrix(0),
                linf=linf, vbk=vbk, t0=t0, M=M, h=h, AgeMax=AgeMax, lbhighs=lbhighs, lbmids=lbmids,
                Mat_a=Mat_a, lwa=lwa, lwb=lwb, 
                Fpen=Fpen, Dpen=Dpen, Dprior=Dprior, 
                RecDev_biasadj=RecDev_biasadj)       
        }
        
        Parameters <- list(log_F_sd=log(SigmaF), log_F_t_input=log(rep(F1,Nyears)),log_q_I=log(qcoef), beta=log(R0), 
            log_sigma_R=log(SigmaR), S50=S50, S95=S95, log_sigma_C=log_sigma_C, log_sigma_I=log_sigma_I, 
            log_CV_L=log_CV_L,Nu_input=rep(0,Nyears))

            Map = list()
            Map[["log_F_sd"]] <- NA
            Map[["log_F_sd"]] <- factor(Map[["log_F_sd"]])
            Map[["S95"]] <- NA
            Map[["S95"]] <- factor(Map[["S95"]])
       #      # Map[["log_sigma_R"]] <- NA
       #      # Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

            if(all(grepl("log_sigma_C", est_sigma))==FALSE){
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])
            }
            if(all(grepl("log_sigma_I", est_sigma))==FALSE){
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
            }
            if(all(grepl("log_CV_L", est_sigma))==FALSE){
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]])
            }
            if(est_sigma==FALSE){
                Map[["log_CV_L"]] <- NA
                Map[["log_CV_L"]] <- factor(Map[["log_CV_L"]]) 
                Map[["log_sigma_C"]] <- NA
                Map[["log_sigma_C"]] <- factor(Map[["log_sigma_C"]])  
                Map[["log_sigma_I"]] <- NA
                Map[["log_sigma_I"]] <- factor(Map[["log_sigma_I"]])
            }

       #  if(RecType==2){
       #      Map[["Nu_input"]] <- 1:length(Parameters[["Nu_input"]])
       #      Map[["Nu_input"]] <- rep(NA, length(Parameters[["Nu_input"]]))
       #      Map[["Nu_input"]] <- factor(Map[["Nu_input"]])

       #      Map[["log_F_t_input"]] <- 1:length(Parameters[["log_F_t_input"]])
       #      Map[["log_F_t_input"]][2:length(Parameters[["log_F_t_input"]])] <- NA
       #      Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
       #  }

       #  if(grepl("Moderate", model)){
       #      Map[["log_sigma_R"]] <- NA
       #      Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])
       #  }

       #  if(grepl("Moderate_Sample", model)){
       #      Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
       #      Map[["log_F_t_input"]][c(1,2,5,6)] <- NA
       #      Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
       #  }

        if(grepl("Poor_Index", model)){
            Map[["log_sigma_R"]] <- NA
            Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

            # Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
            # Map[["log_F_t_input"]][1:10] <- NA
            # Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])

            Map[["beta"]] <- NA
            Map[["beta"]] <- factor(Map[["beta"]])
        }

        if(grepl("Poor_Comp", model)){
            Map[["log_q_I"]] <- NA
            Map[["log_q_I"]] <- factor(Map[["log_q_I"]])

            if(n_lc==1){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][1:10] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
            }
            if(n_lc>1 & n_lc<=10){
                Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
                Map[["log_F_t_input"]][1:10] <- NA
                Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])
            }

            Map[["log_sigma_R"]] <- NA
            Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

            Map[["beta"]] <- NA
            Map[["beta"]] <- factor(Map[["beta"]])
        }

       # if(grepl("Poor_Rel", model)){
       #      Map[["log_q_I"]] <- NA
       #      Map[["log_q_I"]] <- factor(Map[["log_q_I"]])

       #      Map[["beta"]] <- NA
       #      Map[["beta"]] <- factor(Map[["beta"]])

       #      Map[["log_sigma_R"]] <- NA
       #      Map[["log_sigma_R"]] <- factor(Map[["log_sigma_R"]])

       #      # Map[["log_F_t_input"]] = 1:length(Parameters[["log_F_t_input"]])
       #      # Map[["log_F_t_input"]][1] <- NA
       #      # Map[["log_F_t_input"]] <- factor(Map[["log_F_t_input"]])

       #  }

        if(length(Map)==0) Map <- NULL


        if(RecType!=2) Random <- c("Nu_input")
        if(RecType==2) Random <- NULL
        if(REML==TRUE) Random <- union(Random, "log_F_t_input")

    }




	Return <- list("Parameters"=Parameters, "Data"=Data, "Random"=Random, "Map"=Map)
	return(Return)

}