
generateData <- function(modpath, itervec, spatial, Fdynamics, Rdynamics="Constant", LType=1, plotML=FALSE, plotLF_compare=FALSE, plotLF=FALSE){

  ## life history - truth without spatial structure - deterministic across iterations
  lh_nospace <- create_lh_list(lh="Siganus_sutor", selex="asymptotic")

  Fdyn_dir <- file.path(modpath, paste0("F_", Fdynamics))
  dir.create(Fdyn_dir, showWarnings=FALSE)

  for(iter in itervec){

    iterpath <- file.path(Fdyn_dir, iter)
    dir.create(iterpath, showWarnings=FALSE)  

    ## simulated data with no spatial structure in growth
    DataList <- with(lh_nospace, SimData_LB(Nyears=Nyears, AgeMax=AgeMax,
      M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
      SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics,
      R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs,
      lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, Amat=Amat,
      comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears_comp,
      alt_yrs=FALSE, sample=FALSE, nburn=nburn, seed=iter)) 

    # simulated data with spatial structure
    if(spatial==TRUE){

      set.seed(max(itervec)+iter)
      ## spatial process in Linf - varies with each iteration
      spatial_sim <- spatialgrowth_sim(n_i=n_i, Scale=Scale, Sigma2=Sigma2, SD_spatial=SD_spatial, linf=linf, beta_y=beta_y)
      saveRDS(spatial_sim, file.path(iterpath, "spatial_sim.rds"))  

      ## life history - truth with spatial structure - varies with each iteration
      lh_spatial <- lapply(1:nrow(spatial_sim), function(x) create_lh_list(lh="Siganus_sutor", sens="linf", val=spatial_sim[x,"linf_i"], selex="asymptotic")) 

      ## simulated data with spatial structure in growth
      DataList_site <- lapply(1:length(lh_spatial), function(x) with(lh_spatial[[x]], SimData_LB(Nyears=Nyears, AgeMax=AgeMax,
          M=M, F1=F1, h=h, S_a=S_a, qcoef=qcoef, Frate=Frate, Fequil=Fequil, 
          SigmaF=SigmaF, Fdynamics=Fdynamics, Rdynamics=Rdynamics,
          R0=R0, Fmax=Fmax, CVlen=CVlen, mids=mids, highs=highs,
          lows=lows, W_a=W_a, L_a=L_a, Mat_a=Mat_a, Amat=Amat,
          comp_sample=comp_sample, SigmaR=SigmaR, Nyears_comp=Nyears_comp,
          alt_yrs=FALSE, sample=FALSE, nburn=nburn, seed=iter)))  
      SPR_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$SPR)
      RelAbund_site <- sapply(1:length(DataList_site), function(x) DataList_site[[x]]$D_t[length(DataList_site[[x]]$D_t)])
      saveRDS(SPR_site, file.path(iterpath, "SPR_site.rds"))
      saveRDS(RelAbund_site, file.path(iterpath, "RelAbund_site.rds"))

      ## length frequency data at each site
      LF_site <- lapply(1:length(DataList_site), function(x) DataList_site[[x]]$LF)
      ncols_site <- sapply(1:length(LF_site), function(x) ncol(LF_site[[x]]))
      for(i in 1:length(LF_site)){
        ncol <- ncol(LF_site[[i]])
        if(ncol < max(ncols_site)){
          add <- max(ncols_site) - ncol
          LF_site[[i]] <- cbind(LF_site[[i]], matrix(0, nrow=nrow(LF_site[[i]]), ncol=add))
        }
      } 


      ## length frequency by site
      LF_site_array <- array(NA, dim=c(dim(LF_site[[1]]), length(LF_site)))
      ML_t_site <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=length(LF_site))
      for(i in 1:length(LF_site)){
        LF_site_array[,,i] <- LF_site[[i]]
        ML_t_site[,i] <- sapply(1:nrow(LF_site_array[,,i]), function(x) sum(LF_site[[i]][x,]*1:ncol(LF_site[[i]]))/sum(LF_site[[i]][x,]))
      }
      rownames(LF_site_array) <- rownames(ML_t_site) <- (Nyears-Nyears_comp+1):Nyears

      ## length frequency pooled across sites
      LF_pool <- matrix(NA, nrow=nrow(LF_site[[1]]), ncol=max(ncols_site))
      for(i in 1:nrow(LF_site[[1]])){
        for(j in 1:ncol(LF_site[[1]])){
          LF_pool[i,j] <- sum(sapply(1:length(LF_site), function(x) LF_site[[x]][i,j]))
        }
      }
      rownames(LF_pool) <- (Nyears-Nyears_comp+1):Nyears

    if(plotLF==TRUE){
      ## length frequency in the last year at each site
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      
      for(i in 1:length(LF_site)){
        barplot(LF_site[[i]][nrow(LF_site[[i]]),]/sum(LF_site[[i]][nrow(LF_site[[i]]),]), axes=F, xlim=c(0,45), ylim=c(0,0.2))
        mtext(paste0("site ", i), side=3, font=2, line=-3, cex=2)
        if(i %in% 12:15) axis(1, at=pretty(c(0,45)))
        if(i %% 4==0) axis(2, at=pretty(c(0,0.2)))
      }
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)
    }
    if(plotLF_compare==TRUE){
      par(mfrow=c(2,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
      barplot(DataList$LF[nrow(DataList$LF),]/sum(DataList$LF[nrow(DataList$LF),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="tomato")
      mtext(side=3, "no spatial process", font=2, line=-3, cex=2)
      axis(2, at=pretty(c(0,0.2)))      

      barplot(LF_pool[nrow(LF_pool),]/sum(LF_pool[nrow(LF_pool),]), axes=F, xlim=c(0,45), ylim=c(0,0.2), col="black")
      mtext(side=3, "spatial process pooled", font=2, line=-3, cex=2)
      axis(1, at=pretty(c(0,45)))
      axis(2, at=pretty(c(0,0.2)))
      mtext(side=1, "Length bin (1 cm)", outer=TRUE, line=3)
      mtext(side=2, "Proportion", outer=TRUE, line=3)      
    }


      ## mean length over time pooled across sites
      ML_t_pool <- sapply(1:nrow(LF_pool), function(x) sum(LF_pool[x,]*1:ncol(LF_pool))/sum(LF_pool[x,])) 

      if(plotML==TRUE){
      # ## plot mean length at each site
      # # png("SIM_Mean_length_by_site_year.png", res=200, units="in", width=12, height=9)
      par(mfrow=c(4,4), mar=c(0,0,0,0), omi=c(1,1,1,1))
      plot(ML_t_pool, col="red", type="o", pch=17, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
      axis(2, cex=1.2, las=2)
      mtext(side=3, "pooled", font=2, line=-1.5)
      for(i in 1:length(DataList_site)){
        plot(DataList_site[[i]]$ML_t, col="black", pch=19, lwd=2, ylim=c(0,linf), xaxt="n", yaxt="n")
        lines(ML_t_pool, col="red", pch=17, type="o")
        if(i %in% c(12:15)) axis(1, cex=1.2)
        if(i %in% c(4,8,12)) axis(2, cex=1.2, las=2)
        mtext(paste0("site ", i), side=3, font=2, line=-1.5)
      }
      legend("bottomright", legend=c("site-specific", "pooled"), pch=c(19,17), col=c("black", "red"))
      mtext("Year", outer=TRUE, line=3, side=1)
      mtext("Mean length in catch (cm)", outer=TRUE, line=3,  side=2)
      # # dev.off()
    }
      if(LType==1){
        DataList$LF <- LF_pool
        DataList$ML_t <- ML_t_pool  
      }
      if(LType==0){
        DataList$LF <- LF_site_array
        DataList$ML_t <- ML_t_site
      }


      rm(spatial_sim)
      rm(lh_spatial)
    }

    if(ncol(DataList$LF) > 50) stop("More than 50 length bins - must resize")
    DataList_out <- DataList
    if(ncol(DataList$LF) < 50){
      add <- 50 - ncol(DataList$LF)
      DataList_out$LF <- cbind(DataList$LF, matrix(0, nrow=nrow(DataList$LF), ncol=add))
    }

      saveRDS(DataList_out, file.path(iterpath, "True.rds"))
      rm(DataList)
      rm(DataList_out)
      rm(iterpath)

}

  return(paste0(length(itervec), " iterates of data generated in ", modpath))

}