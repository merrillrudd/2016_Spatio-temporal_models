calc_Fref <- function(Mat_a, W_a, M, S_a, R0, F, ref){

		Na0 <- Naf <- rep(NA, length(W_a))
		Na0[1] <- Naf[1] <- 1
		for(a in 2:length(W_a)){
			if(a<length(W_a)){
				Na0[a] <- Na0[a-1]*exp(-M)
				Naf[a] <- Naf[a-1]*exp(-M-S_a[a-1]*F)
			}
			if(a==length(W_a)){
				Na0[a] <- (Na0[a-1]*exp(-M))/(1-exp(-M))
				Naf[a] <- (Naf[a-1]*exp(-M-F*S_a[a-1]))/(1-exp(-M-F*S_a[a-1]))
			}
		}

		SB0 <- sum(Na0[-1]*Mat_a[-1]*W_a[-1])
		SBf <- sum(Naf[-1]*Mat_a[-1]*W_a[-1])

		SPR <- SBf/SB0
		if(ref==FALSE) return(SPR)
		
		if(ref!=FALSE){
			diff <- ref - depl
			return(diff)
		}
}