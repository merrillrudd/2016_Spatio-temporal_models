rm(list=ls())

data_dir <- "C:\\Git_Projects\\2016_Spatio-temporal_models\\Project"
# data_dir <- "F:\\Merrill\\Git_Projects\\2016_Spatio-temporal_models\\Project"
setwd(data_dir)
source("R_functions\\functions.R")

require(geoR)
require(RandomFields)
require(ape)
require(RColorBrewer)
require(spdep)
require(reshape2)
require(TMB)

### QUESTION
### Is there bias in estimates of stock status when there is spatial variation in growth, but it is ignored in the assessment?

### PERFORMANCE MEASURE
#### the ability to estimate the mean SPR across sites with spatial variation
#### or how many sites is SPR over-estimated

### METHOD STRUCTURE
## OM1 = no spatial differences in growth, constant recruitment with process error, constant/endogenous fishing mortality with process error
## OM1 = spatial differences in growth, constant recruitment with process error, constant/endogenous fishing mortality with process error
## EM1 = growth is the same across space

### EXECUTABLE
run_exe <- file.path(data_dir, "inst", "executables")

Version <- "lb_statespace_v5"
setwd(run_exe)
dyn.unload( paste0(Version,".dll") )
file.remove( paste(Version,c(".dll",".o"),sep="") )
compile(paste0(Version, ".cpp"))

### RESULTS DIRECTORIES
pool_space_dir <- file.path(data_dir, "pool_space")
dir.create(pool_space_dir, showWarnings=FALSE) 

unlink(pool_space_dir, TRUE)
dir.create(pool_space_dir, showWarnings=FALSE) 

nospace_dir <- file.path(pool_space_dir, "nospace_dir")
dir.create(nospace_dir, showWarnings=FALSE)

space_dir <- file.path(pool_space_dir, "space_dir")
dir.create(space_dir, showWarnings=FALSE)

### SETTINGS - GENERAL
dat_scenario <- "Poor_Comp_LC"
Nyears <- 20
comp_sample <- 100 
SigmaR <- 0.6
seed <- 123
Nyears_comp <- 10
simdir <- "pool_space"
sens_param <-FALSE
sens_val <- FALSE
RecType <- 0
RecTypeDesc <- "shared"
OMselex <- "asymptotic"
EMselex <- "asymptotic"
nburn <- 20
est_sigma <- FALSE
rerun=FALSE
runLBSPR=FALSE
REML=FALSE

### SETTINGS - SPATIAL
### simulation of 1D spatial process on growth - from HW5
n_i <- 15  ## number of sites
Scale <- 2 ## scale for Gaussian random field
Sigma2 <- 1 ## variance for random effect
SD_spatial <- exp(0.1) ## variance for Gaussian random field
linf <- 36.2  ## Linf_0
beta_y <- 0.02 ## spatial trend

setwd(data_dir)
source("R_functions\\functions.R")

### MODEL PATHS
modpath_nospace <- file.path(data_dir, "pool_space", "nospace_dir")
modpath_space <- file.path(data_dir, "pool_space", "space_dir")
modpath_vec <- c(modpath_nospace, modpath_space)
mod_names <- c("No spatial process", "Spatial process pooled")

### DATA INPUTS
## life history - truth without spatial structure - deterministic across iterations
lh_nospace <- create_lh_list(lh="Siganus_sutor", selex="asymptotic")
## data inputs - same across scenarios, assumes no spatial differences in growth
dat_input <- create_inputs(param=FALSE, val=FALSE, lh_dat=lh_nospace)

### RUN MODEL
## number of simulations
niter <- 1

start <- Sys.time()

## loop over fishing mortality patterns (built into OM)
Fdyn_vec <- c("Constant", "Endogenous", "Ramp")

## generate data - no spatial structure in true population
sapply(1:length(Fdyn_vec), function(x) generateData(modpath=modpath_nospace, itervec=1:niter, spatial=FALSE, Fdynamics=Fdyn_vec[x]))

## generate data - spatial structure of linf in true population
sapply(1:length(Fdyn_vec), function(x) generateData(modpath=modpath_space, itervec=1:niter, spatial=TRUE, Fdynamics=Fdyn_vec[x]))

## run model - no spatial structure in true population
sapply(1:length(Fdyn_vec), function(x) runModel(modpath=modpath_nospace, itervec=1:niter, data_input=dat_input, Fdynamics=Fdyn_vec[x]))

## run model - spatial structure of linf in true population
sapply(1:length(Fdyn_vec), function(x) runModel(modpath=modpath_space, itervec=1:niter, data_input=dat_input, Fdynamics=Fdyn_vec[x]))

end <- Sys.time() - start

### FIGURES
## results - relative error SPR
par(mfrow=c(3,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
res_constant <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Constant", mod_names=mod_names)
res_endogenous <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Endogenous", mod_names=mod_names)
res_ramp <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Ramp", mod_names=mod_names)
axis(1, at=c(1,2), labels=mod_names, cex.axis=2)
mtext("Relative Error", side=2, line=2, cex=1.3, outer=TRUE)

## mean length across sites - generated data
plotML <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotML=TRUE)
## length composition across sites
plotLF <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotLF=TRUE)
## compare length comp with no spatial variation to pooled length comp with spatial variation 
plotLFcompare <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotLF_compare=TRUE)

## asymptotic length across sites
linf_sites <- y_sites <- matrix(NA, nrow=n_i, ncol=64)
for(i in 1:64){
	sim <- readRDS(file.path(modpath_space, "F_Constant", i, "spatial_sim.rds"))
	linf_sites[,i] <- sim$linf_i
	y_sites[,i] <- sim$y_i
}
dim <- 5
par(mfrow=c(dim,dim), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:(dim^2)){
	plot(y_sites[,i], linf_sites[,i],axes=F, ann=F, ylim=c(min(linf_sites), max(linf_sites)), xlim=c(-4,-1),pch=19)
	abline(h=linf, col="red", lty=2)
	if(i==1 | i %% dim == 1) axis(2, at=pretty(c(min(linf_sites), max(linf_sites))), las=2)
	if(i %in% (dim*(dim-1)+1):(dim^2)) axis(1, at=seq(-3.5,-1.5,by=0.5))
	text(x=-3.9, y=max(linf_sites)*0.97, letters[i], font=2, cex=1.5)
	box(col="gray")
}
mtext(side=1, "Latitude", outer=TRUE, line=3)
mtext(side=2, "Asymptotic length (cm)", outer=TRUE, line=3)

## compare model fits
plotFIT(compare_quant=c( "ML", "R", "F", "D"), compare_type="base_values", scenario_name="No spatial structure", modpath=file.path(modpath_nospace, "F_Constant"), iter=1)
dev.new()
plotFIT(compare_quant=c( "ML", "R", "F", "D"), compare_type="base_values", scenario_name="Spatial structure", modpath=file.path(modpath_space, "F_Constant"), iter=1)


