create_lh_list <- function(lh, sens=FALSE, val=FALSE, selex){

	if(lh==1){
		vbk <- 0.2
		linf <- 37.3
		t0 <- -0.2
		binwidth <- 1
		CVlen <- 0.1
		R0 <- 1
		lwa <- 0.0214
		lwb <- 2.94
		h <- 0.7
		Fequil <- 0.25
		Frate <- 0.2
		F1 <- 0.2
		qcoef <- 1e-2
		Fmax <- 1
		SigmaF <- 0.3
		SigmaC <- 0.2
		SigmaI <- 0.2
		mids <- seq((binwidth/2), linf*1.2, by=binwidth)
		highs <- mids + (binwidth/2)
		lows <- mids - (binwidth)/2
		L50 <- floor((2/3)*linf)
		M <- 1.84*vbk

		S50 <- 1
		S95 <- 2
		SL50 <- round(linf*(1-exp(-vbk*(S50-t0))))
		SL95 <- round(linf*(1-exp(-vbk*(S95-t0))))
	}

	if(lh=="Siganus_sutor"){
		vbk <- 0.87
		if(sens=="vbk") vbk <- val
		linf <- 36.2
		if(sens=="linf") linf <- val
		t0 <- -0.2
		binwidth <- 1
		CVlen <- 0.1
		if(sens=="CVlen") CVlen <- val
		R0 <- 1
		lwa <- 0.05970
		lwb <- 2.754
		h <- 1
		Fequil <- 0.25
		Frate <- 0.2
		F1 <- 0.2
		qcoef <- 1e-2
		Fmax <- 1
		SigmaF <- 0.3
		SigmaC <- 0.2
		SigmaI <- 0.2
		mids <- seq((binwidth/2), 50, by=binwidth)
		highs <- mids + (binwidth/2)
		lows <- mids - (binwidth)/2
		
		M <- 1.6*vbk
		L50 <- 20.2
		SL50 <- 11.3
		S50 <- ceiling(t0-log(1-(SL50/linf))/vbk)
		S95 <- S50+1
		SL95 <- ceiling(linf*(1-exp(-vbk*(S95-t0))))
	}


	AgeMax <- max(round(-log(0.01)/M), 6)
	ages <- 0:AgeMax
	Amat <- round(t0-log(1-(L50/linf))/vbk)
	A95 <- Amat+1
	L95 <- round(linf*(1-exp(-vbk*(A95-t0))))

	L_a <- linf*(1-exp(-vbk*(ages - t0)))
	W_a <- lwa*L_a^lwb
	# Amat <- log( 3*(linf - 10)/linf )/vbk
	Mat_a <- 1 / (1 + exp(Amat - ages))
	if(selex=="asymptotic"){
		S_a <- 1/(1+exp(-log(19)*(ages-S50)/(S95-S50))) # Selectivity at age
		Syoung <- NA
		Sold <- NA
	}
	if(selex=="dome"){
		S_a_calc <- rep(NA, length(ages))
		Syoung <- 0.8
		Sold <- 8
		A <- sqrt(2/pi)/(Syoung + Sold)
		for(a in 1:length(ages)){
			if(a <= S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Syoung^2))
			if(a > S95) S_a_calc[a] <- A*exp(-((ages[a] - S95)^2)/(2*Sold^2))
		}
		S_a <- S_a_calc/max(S_a_calc)
	}
	


	Outs <- NULL
	Outs$vbk <- vbk
	Outs$linf <- linf
	Outs$t0 <- t0
	Outs$binwidth <- binwidth
	Outs$CVlen <- CVlen
	Outs$SigmaC <- SigmaC
	Outs$SigmaI <- SigmaI
	Outs$R0 <- R0
	Outs$lwa <- lwa
	Outs$lwb <- lwb
	Outs$S50 <- S50
	Outs$S95 <- S95
	Outs$SL50 <- SL50
	Outs$SL95 <- SL95
	Outs$Syoung <- Syoung
	Outs$Sold <- Sold
	Outs$h <- h
	Outs$Fequil <- Fequil
	Outs$Frate <- Frate
	Outs$F1 <- F1
	Outs$qcoef <- qcoef
	Outs$Fmax <- Fmax
	Outs$SigmaF <- SigmaF
	Outs$M <- M
	Outs$AgeMax <- AgeMax
	Outs$mids <- mids
	Outs$highs <- highs
	Outs$lows <- lows
	Outs$S_a <- S_a
	Outs$L_a <- L_a
	Outs$W_a <- W_a
	Outs$Amat <- Amat
	Outs$L50 <- L50
	Outs$L95 <- L95
	Outs$Mat_a <- Mat_a

	return(Outs)

}