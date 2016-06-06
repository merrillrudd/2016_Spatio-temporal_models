compare_LBSPR <- function(lh_num, Fdynamics, Rdynamics, itervec, simdir,
	SigmaR, test_dir, dat_scenario, Ncomp, sens_param, sens_label, RecTypeDesc, 
	OMselex, EMselex, year="terminal"){

	sprF <- ssaF <- sprD <- ssaD <- trueF <- trueD <- sprSPR <- ssaSPR <- trueSPR <- rep(NA, length(itervec))

	for(iter in itervec){
		iterpath <- setPaths(lh=lh_num, Fdyn=Fdynamics, Rdyn=Rdynamics, iter=iter, SimDir=simdir,
        SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, iterpath=TRUE, 
        Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc,
        OMselex=OMselex, EMselex=EMselex)

		if(simdir=="base") lbspr_res <- readRDS(file.path(iterpath, "LBSPR_results.rds"))
		if(simdir=="explore_comp") lbspr_res <- readRDS(file.path(iterpath, paste0("LBSPR_results_", (20-year +1), ".rds")))
		lbssa_res <- readRDS(file.path(iterpath, "Report.rds"))
		truth <- readRDS(file.path(iterpath, "True.rds"))

		sprF[which(itervec==iter)] <- lbspr_res$F
		sprD[which(itervec==iter)] <- lbspr_res$depl
		sprSPR[which(itervec==iter)] <- lbspr_res$SPR

		if(year=="terminal") ssaF[which(itervec==iter)] <- lbssa_res$F_t[length(lbssa_res$F_t)]
		if(year!="terminal") ssaF[which(itervec==iter)] <- lbssa_res$F_t[year]

		if(year=="terminal") ssaD[which(itervec==iter)] <- lbssa_res$Depl[length(lbssa_res$Depl)]
		if(year!="terminal") ssaD[which(itervec==iter)] <- lbssa_res$Depl[year]

		if(year=="terminal") ssaSPR[which(itervec==iter)] <- calc_ref(lbssa_res$Mat_a, lbssa_res$W_a, lbssa_res$M, lbssa_res$S_a, lbssa_res$R0, lbssa_res$F_t[length(lbssa_res$F_t)], ref=FALSE)
		if(year!="terminal") ssaSPR[which(itervec==iter)] <- calc_ref(lbssa_res$Mat_a, lbssa_res$W_a, lbssa_res$M, lbssa_res$S_a, lbssa_res$R0, lbssa_res$F_t[year], ref=FALSE)

		if(year=="terminal") trueF[which(itervec==iter)] <- truth$F_t[length(truth$F_t)]
		if(year!="terminal") trueF[which(itervec==iter)] <- truth$F_t[year]

		if(year=="terminal") trueD[which(itervec==iter)] <- truth$D_t[length(truth$D_t)]
		if(year!="terminal") trueD[which(itervec==iter)] <- truth$D_t[year]

		if(year=="terminal") trueSPR[which(itervec==iter)] <- calc_ref(truth$Mat_a, truth$W_a, truth$M, truth$S_a, truth$R0, truth$F_t[length(truth$F_t)], ref=FALSE)
		if(year!="terminal") trueSPR[which(itervec==iter)] <- calc_ref(truth$Mat_a, truth$W_a, truth$M, truth$S_a, truth$R0, truth$F_t[year], ref=FALSE)


	}

	

	Outs <- NULL
	Outs$sprF <- sprF
	Outs$sprD <- sprD
	Outs$ssaF <- ssaF
	Outs$ssaD <- ssaD
	Outs$trueF <- trueF
	Outs$trueD <- trueD
	Outs$sprSPR <- sprSPR 
	Outs$ssaSPR <- ssaSPR
	Outs$trueSPR <- trueSPR

	return(Outs)

}