LBSPR_sim <- function(lh_list, setSPR){

	Linf <- lh_list$linf
	vbk <- lh_list$vbk
	t0 <- lh_list$t0

	L50 <- lh_list$Lmat

	S50 <- lh_list$S50
	S95 <- lh_list$S95
	SL50 <- round(Linf*(1-exp(-vbk*(S50-t0))))
	SL95 <- round(Linf*(1-exp(-vbk*(S95-t0))))

	M <- lh_list$M

	lens <- seq(from=0, to=max(L50, Linf), length.out=100)
	Matprob <-	1.0/(1+exp((lens-L50)))
	Selprob <-	1.0/(1+exp(-log(19)*(lens-SL50)/(SL95-SL50)))
	
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
    SizeBins$Linc <- lh_list$binwidth
    SizeBins$ToSize <- Linf*1.25

	# Calc SPR v F/M 
    FMVec <- seq(from=0, to=5, length.out=100)
    run <- sapply(1:length(FMVec), function (xx) {
      Fleet$FM <- FMVec[xx]
      EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)
    })
	SaveSPR <- sapply(1:length(FMVec), function(X) run[,X]$SPR)
	SaveYield1 <- sapply(1:length(FMVec), function(X) run[,X]$Yield)
    SaveYield <- SaveYield1/max(SaveYield1)
    FTarg <- FMVec[max(which(SaveSPR >= setSPR))]
  	
	# getFM <- optimise(getFMfun, interval=c(0, 100), setSPR=setSPR, Stock=Stock, Fleet=Fleet, SizeBins=SizeBins, FitControl=NULL)
	# print(getFM)
	Fleet$FM <- FTarg
    runMod <- EqSimMod_LB(Stock, Fleet, SizeBins, FitControl=NULL)
	currYield <- runMod$Yield/max(SaveYield1)
	
}