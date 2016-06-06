plot_LBSPR <- function(itervec, simdir, SigmaR, test_dir, dat_scenario="Poor_Comp_LC", Ncomp=NULL, sens_param=FALSE, sens_label=NULL, RecTypeDesc, OMselex, EMselex, year="terminal"){

if(year=="terminal"){
	par(mfrow=c(3,3), mar=c(0,0,0,0), omi=c(1,1,1,1))
	compare <- compare_LBSPR(lh_num=1, Fdynamics="Endogenous", Rdynamics="Constant", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", las=2)
	abline(h=0, lty=2, lwd=2)
	mtext(side=3, "F Endogenous", font=2, line=2)
	mtext(side=2, "R Constant", font=2, line=3)

	compare <- compare_LBSPR(lh_num=1, Fdynamics="Constant", Rdynamics="Constant", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")
	
	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")	
	abline(h=0, lty=2, lwd=2)
	mtext(side=3, "F Constant", font=2, line=2)

	compare <- compare_LBSPR(lh_num=1, Fdynamics="Ramp", Rdynamics="Constant", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")	
	abline(h=0, lty=2, lwd=2)
	mtext(side=3, "F Ramp", font=2, line=2)

	compare <- compare_LBSPR(lh_num=1, Fdynamics="Endogenous", Rdynamics="Pulsed", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", las=2)	
	abline(h=0, lty=2, lwd=2)
	mtext(side=2, "R Pulsed", font=2, line=3)

	compare <- compare_LBSPR(lh_num=1, Fdynamics="Constant", Rdynamics="Pulsed", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")	
	abline(h=0, lty=2, lwd=2)

	compare <- compare_LBSPR(lh_num=1, Fdynamics="Ramp", Rdynamics="Pulsed", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")
	abline(h=0, lty=2, lwd=2)


	compare <- compare_LBSPR(lh_num=1, Fdynamics="Endogenous", Rdynamics="BH", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", las=2)
	abline(h=0, lty=2, lwd=2)
	mtext(side=2, "R BH", font=2, line=3)
	axis(1, at=c(1,2), labels=c("LB-SPR", "LB-SSA"))


	compare <- compare_LBSPR(lh_num=1, Fdynamics="Constant", Rdynamics="BH", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")
	abline(h=0, lty=2, lwd=2)
	axis(1, at=c(1,2), labels=c("LB-SPR", "LB-SSA"))


	compare <- compare_LBSPR(lh_num=1, Fdynamics="Ramp", Rdynamics="BH", itervec=itervec, simdir=simdir,  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=Ncomp, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year="terminal")

	relerr <- list()
    relerr$sprF <- (compare$sprF - compare$trueF)/compare$trueF
    relerr$sprD <- (compare$sprD - compare$trueD)/compare$trueD
    relerr$ssaF <- (compare$ssaF - compare$trueF)/compare$trueF
    relerr$ssaD <- (compare$ssaD - compare$trueD)/compare$trueD

    relerr$sprSPR <- (compare$sprSPR - compare$trueSPR)/compare$trueSPR
    relerr$ssaSPR <- (compare$ssaSPR - compare$trueSPR)/compare$trueSPR

	# boxplot(relerr$sprF, relerr$ssaF, relerr$sprD, relerr$ssaD, col=c("violet", "goldenrod", "violet", "goldenrod"), ylim=c(-1,10))
	boxplot(relerr$sprSPR, relerr$ssaSPR, ylim=c(-1, 5), col=c("violet", "goldenrod"), xaxt="n", yaxt="n")
	abline(h=0, lty=2, lwd=2)
	axis(1, at=c(1,2), labels=c("LB-SPR", "LB-SSA"))
}

if(year!="terminal"){

par(mfrow=c(3,3))
	iterpath <- setPaths(lh=1, Fdyn="Endogenous", Rdyn="Constant", iter=1, SimDir=simdir,
        SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, iterpath=TRUE, 
        Ncomp=5, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc,
        OMselex=OMselex, EMselex=EMselex)

		truth <- readRDS(file.path(iterpath, "True.rds"))
		plot(truth$D_t, lwd=2, type="l")

		compare <- compare_LBSPR(lh_num=1, Fdynamics="Endogenous", Rdynamics="Constant", itervec=itervec, simdir="explore_comp",  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=5, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year=20)

		points(x=20, compare$sprD[1], cex=1.2, pch=18, col="violet")
		points(x=20, compare$ssaD[1], cex=1.2, pch=16, col="goldenrod")

		compare <- compare_LBSPR(lh_num=1, Fdynamics="Endogenous", Rdynamics="Constant", itervec=itervec, simdir="explore_comp",  SigmaR=SigmaR, test_dir=test_dir, dat_scenario=dat_scenario, Ncomp=5, sens_param=sens_param, sens_label=sens_label, RecTypeDesc=RecTypeDesc, OMselex=OMselex, EMselex=EMselex, year=19)

		points(x=19, compare$sprD[1], cex=1.2, pch=18, col="violet")
		points(x=19, compare$ssaD[1], cex=1.2, pch=16, col="goldenrod")

}

}