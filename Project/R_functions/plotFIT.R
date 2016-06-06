plotFIT <- function(compare_quant, compare_type, ylab=FALSE, scenario_name, modpath, iter){
        FUN = function(InputMat, log=TRUE){
          if(log==TRUE) return(c( exp(InputMat[,1]-1.96*InputMat[,2]), rev(exp(InputMat[,1]+1.96*InputMat[,2]))))
          if(log==FALSE) return(c( InputMat[,1]-1.96*InputMat[,2], rev(InputMat[,1]+1.96*InputMat[,2])))
        } 
if(compare_type=="base_values"){
    ncol <- length(compare_quant)/2
    nrow <- 2
}
if(compare_type=="base_scenarios"){
    ncol <- length(modpath)
    nrow <- length(compare_quant)
}
par(mfrow=c(nrow, ncol), mar=c(0,3,0,0), omi=c(1,1,1,1))

for(dd in 1:length(modpath)){
  True <- readRDS(file.path(modpath[dd], iter, "True.rds"))
  Report <- readRDS(file.path(modpath[dd], iter, "Report.rds"))
  Sdreport <- readRDS(file.path(modpath[dd], iter, "Sdreport.rds"))

  for(qq in 1:length(compare_quant)){
          if("B"==compare_quant[qq]){
            ## Biomass
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$SB_t, "Est"=Report$SB_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lSB_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Biomass", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Biomass", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

          }
          if("N"==compare_quant[qq]){
            ## Abundance
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$N_t, "Est"=Report$N_t_hat)
            ymax <- 4
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lN_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Abundance", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Abundance", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

          }
          if("ML"==compare_quant[qq]){
            ## Average Length
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$L_t, "Est"=Report$L_t_hat)
            ymax <- 40
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="L_t_hat"),], log=FALSE), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Mean Length", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Mean Length", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

          }
          if("R"==compare_quant[qq]){       
            ## Recruitment
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$R_t, "Est"=Report$R_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE) if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lR_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Recruitment", side=3, line=-2, font=2)
            if(ylab==TRUE) mtext("Recruitment", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

          }       
          if("F"==compare_quant[qq]){
            ## Fishing Mortality
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$F_t, "Est"=Report$F_t)
            ymax <- 2
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lF_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Fishing Mortality", side=3, line=-2, font=2)
            if(compare_type=="base_values" | "F"==compare_quant[length(compare_quant)]) axis(1, at=seq(5,20, by=5))   
            if(ylab==TRUE) mtext("Fishing Mortality", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

          }
          if("D"==compare_quant[qq]){     
            ## Depletion
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=True$D_t, "Est"=Report$Depl)
            ymax <- 3
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lD_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Depletion", side=3, line=-2, font=2)  
            if(compare_type=="base_values"| "D"==compare_quant[length(compare_quant)]) axis(1, at=seq(5,20, by=5))   
            if(ylab==TRUE) mtext("Depletion", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)          
          }
          if("C"==compare_quant[qq]){      
            ## Catch
            tcatch <- True$C_t
            pcatch <- rep(NA, length(True$SB_t))
              names(pcatch) <- 1:length(True$SB_t)
            pcatch[which(names(pcatch) %in% names(True$C_t))] <- tcatch
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=pcatch, "Est"=Report$C_t_hat)
            ymax <- 5
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lC_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Catch", side=3, line=-2, font=2)   
            if(compare_type=="base_values"| "C"==compare_quant[length(compare_quant)]) axis(1, at=seq(5,20, by=5))   
            if(ylab==TRUE) mtext("Catch", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)        
          }
          if("I"==compare_quant[qq]){     
            ## Index
            tindex <- True$I_t
            pindex <- rep(NA, length(True$SB_t))
              names(pindex) <- 1:length(True$SB_t)
            pindex[which(names(pindex) %in% names(True$I_t))] <- tindex
            Mat <- cbind("Year"=1:length(True$SB_t), "True"=pindex, "Est"=Report$I_t_hat)
            ymax <- 3
            matplot(y=Mat[,c("True", "Est")], x=Mat[,c("Year")], type="l", col=c("black", "red"), lty="solid", ylim=c(0, ymax), lwd=2, xaxt="n")
            if(all(is.na(Sdreport))==FALSE)if( !("condition" %in% names(attributes(Sdreport)))) polygon( y=FUN(summary(Sdreport)[which(rownames(summary(Sdreport))=="lI_t"),]), x=c(Mat[,c("Year")], rev(Mat[,c("Year")])), col=rgb(1,0,0,alpha=0.2), border=NA)
            if(dd==1) mtext("Index", side=3, line=-2, font=2)
            if(compare_type=="base_values"| "I"==compare_quant[length(compare_quant)])  axis(1, at=seq(5,20, by=5))           
            if(ylab==TRUE) mtext("Index", side=2, line=3, font=2)
            if(qq==1) mtext(scenario_name, side=3, line=1)

        }
      }
}
      

}