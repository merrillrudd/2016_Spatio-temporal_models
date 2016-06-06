LBSPR_assessment <- function(lh_list, lendat, setSPR){

	Linf <- lh_list$linf
	vbk <- lh_list$vbk
	t0 <- lh_list$t0

	L50 <- lh_list$L50
	L95 <- lh_list$L95

	SL50 <- lh_list$SL50
	SL95 <- lh_list$SL95

	M <- lh_list$M

	binwidth <- lh_list$binwidth

	MaxLen <- round(max(1.25*Linf, 1.1*max(lendat)), 0)
	LenBins <- seq(from=0, to=MaxLen, by=binwidth)
	LenMids <- seq(from=0.5*binwidth, by=binwidth, length.out=length(LenBins)-1)
	LenDat <- as.vector(table(cut(lendat, LenBins)))

	## same as LBSPR_sim
	Stock <- NULL
	Stock$NGTG <- 41
	Stock$GTGLinfBy <- NA
    Stock$Linf <- Linf
    Stock$CVLinf <- 0.1 # NEED TO ADD THIS TO INPUT VARIABLES
    Stock$MaxSD <- 2 
    Stock$MK  <- M/vbk
    Stock$L50 <- L50
    Stock$Walpha <- lh_list$lwa
    Stock$Wbeta <- lh_list$lwb
    Stock$FecB  <- 3
    Stock$Steepness <- lh_list$h
    Stock$Mpow <- 0
    Stock$R0  <- lh_list$R0
    Stock$CentLinf <- NULL
    Stock$CentMpar <- NULL
    Stock$CentKpar <- NULL
    Stock$Mslope <- 0 

    ## same as LBSPR_sim
	Fleet <- NULL
    Fleet$SL50 <- SL50
    Fleet$SL95 <- SL95
    Fleet$MLLKnife <- NA
    Fleet$FM <- 0

	SizeBins <- NULL
    SizeBins$Linc <- binwidth
    SizeBins$ToSize <- max(LenBins)

	SL50Start <- LenMids[which.max(LenDat)]
	DeltaStart <- 0.1*SL50Start
	FMStart <- 1
	Starts <- log(c(SL50Start/Linf, DeltaStart/Linf, FMStart))
	Lower <- log(c(0.1, 0.1, 0.001))
	Upper <- log(c(0.9, 0.9, 20))

	Opt <- nlminb(Starts, OptRoutine, LenDat=LenDat, Stock=Stock, 
		SizeBins=SizeBins, lower=Lower, upper=Upper)
	N <- sum(LenDat)

	Fleet <- NULL
	Fleet$SL50 <- exp(Opt$par)[1] * Stock$Linf
	Fleet$SL95 <- Fleet$SL50 + exp(Opt$par)[2] * Stock$Linf
	Fleet$MLLKnife <- NA
	Fleet$FM <- exp(Opt$par)[3]
	runMod <- EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)

	# Calc SPR v F/M 
    FMVec <- seq(from=0, to=5, length.out=100)
    run <- sapply(1:length(FMVec), function (xx) {
      Fleet$FM <- FMVec[xx]
      EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)
    })
	SaveSPR <- sapply(1:length(FMVec), function(X) run[,X]$SPR)
	SaveYield1 <- sapply(1:length(FMVec), function(X) run[,X]$Yield)
    SaveYield <- SaveYield1/max(SaveYield1)
	currYield <- runMod$Yield/max(SaveYield1)
    FTarg <- FMVec[max(which(SaveSPR >= setSPR))]
	  
	lens <- seq(from=0, to=max(L50, Linf), length.out=100)
	Matprob <-	1.0/(1+exp(-log(19)*(lens-L50)/(L95-L50)))
	Selprob <-	1.0/(1+exp(-log(19)*(lens-SL50)/(SL95-SL50)))

	# MCex <- 1.3
	# par(mfrow=c(2,1))
	#   plot(c(0,5), c(0,1), type="n", bty="n", las=1, xlab="", ylab="")
	#   lines(FMVec, SaveSPR, lwd=3)
	#   mtext(side=1, line=2.5, expression(italic(F/M)), cex=MCex)
	#   mtext(side=2, line=2.5, "SPR", cex=MCex)
	#   points(Fleet$FM, runMod$SPR, pch=19, cex=3, xpd=NA)
	 
	#   plot(c(0,5), c(0,1), type="n", bty="n", las=1, xlab="", ylab="")
	#   lines(FMVec, SaveYield, lwd=3)
	#   mtext(side=1, line=2.5, expression(italic(F/M)), cex=MCex)
	#   mtext(side=2, line=2.5, "Relative Yield", cex=MCex)
	#   points(Fleet$FM, currYield, pch=19, cex=3, xpd=NA)

	 Outs <- NULL
	 Outs$SPR <- runMod$SPR
	 Outs$depl <- runMod$depl
	 Outs$F <- (Fleet$FM)*M
	 return(Outs)
}