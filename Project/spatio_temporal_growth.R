rm(list=ls())

# data_dir <- "C:\\Git_Projects\\2016_Spatio-temporal_models\\Project"
data_dir <- "F:\\Merrill\\Git_Projects\\2016_Spatio-temporal_models\\Project"
setwd(data_dir)
source("R_functions\\functions.R")

require(geoR)
require(RandomFields)
require(ape)
require(RColorBrewer)
require(spdep)
require(reshape2)
require(TMB)

############ Part 1 ###############
### QUESTION
### Is there bias in estimates of stock status when there is spatial variation in growth, but it is ignored in the assessment?

### PERFORMANCE MEASURE
#### the ability to estimate the mean SPR across sites with spatial variation
#### or how many sites is SPR over-estimated

### METHOD STRUCTURE
## OM1 = no spatial differences in growth, constant recruitment with process error, constant/endogenous fishing mortality with process error
## OM1 = spatial differences in growth, constant recruitment with process error, constant/endogenous fishing mortality with process error
## EM1 = growth is the same across space

## locate and compile executable
run_exe <- file.path(data_dir, "inst", "executables")

Version <- "lb_statespace_v5"
setwd(run_exe)
dyn.unload( paste0(Version,".dll") )
file.remove( paste(Version,c(".dll",".o"),sep="") )
compile(paste0(Version, ".cpp"))

## directories for results
pool_space_dir <- file.path(data_dir, "pool_space")
dir.create(pool_space_dir, showWarnings=FALSE) 

unlink(pool_space_dir, TRUE)
dir.create(pool_space_dir, showWarnings=FALSE) 

nospace_dir <- file.path(pool_space_dir, "nospace_dir")
dir.create(nospace_dir, showWarnings=FALSE)

space_dir <- file.path(pool_space_dir, "space_dir")
dir.create(space_dir, showWarnings=FALSE)

## settings - general
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

## settings - space
### simulation of 1D spatial process on growth - from HW5
n_i <- 15 ##  number of sites
Scale <- 2
Sigma2 <- 0.1
SD_spatial <- 0.1
linf <- 36.2 
beta_y <- 0.05

setwd(data_dir)
source("R_functions\\functions.R")

##### demonstrate example of spatial process
# spatial process in Linf - incorporates anabolism and catabolism
spatial_sim <- spatialgrowth_sim(n_i=n_i, Scale=Scale, Sigma2=Sigma2, SD_spatial=SD_spatial, linf=linf, beta_y=beta_y)

## plot linf across latitude
par(mfrow=c(1,1), mar=c(4,4,4,4))
plot(spatial_sim$y_i, spatial_sim$linf_i, xlab="Latitude", ylab="Linf value")

## model path
modpath_nospace <- file.path(data_dir, "pool_space", "nospace_dir")
modpath_space <- file.path(data_dir, "pool_space", "space_dir")
modpath_vec <- c(modpath_nospace, modpath_space)
mod_names <- c("No spatial process", "Spatial process pooled")

## life history - truth without spatial structure - deterministic across iterations
lh_nospace <- create_lh_list(lh="Siganus_sutor", selex="asymptotic")
## data inputs - same across scenarios, assumes no spatial differences in growth
dat_input <- create_inputs(param=FALSE, val=FALSE, lh_dat=lh_nospace)

niter <- 1000

start <- Sys.time()

Fdyn_vec <- c("Constant", "Endogenous", "Ramp")

plotML <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotML=TRUE)
plotLF <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotLF=TRUE)
plotLFcompare <- generateData(modpath=modpath_space, itervec=1, spatial=TRUE, Fdynamics="Constant", LType=0, plotLF_compare=TRUE)

## generate data - no spatial structure in true population
sapply(1:length(Fdyn_vec), function(x) generateData(modpath=modpath_nospace, itervec=1:niter, spatial=FALSE, Fdynamics=Fdyn_vec[x]))

## generate data - spatial structure of linf in true population
sapply(1:length(Fdyn_vec), function(x) generateData(modpath=modpath_space, itervec=1:niter, spatial=TRUE, Fdynamics=Fdyn_vec[x], LType=0))

## run model - no spatial structure in true population
sapply(1:length(Fdyn_vec), function(x) runModel(modpath=modpath_nospace, itervec=1:niter, data_input=dat_input, Fdynamics=Fdyn_vec[x]))

## run model - spatial structure of linf in true population
sapply(1:length(Fdyn_vec), function(x) runModel(modpath=modpath_space, itervec=1:niter, data_input=dat_input, Fdynamics=Fdyn_vec[x]))


end <- Sys.time() - start

par(mfrow=c(3,1), mar=c(0,0,0,0), omi=c(1,1,1,1))
res_constant <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Constant", mod_names=mod_names)
res_endogenous <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Endogenous", mod_names=mod_names)
res_ramp <- SPRerror(modpath_vec=modpath_vec, niter=niter, Fdyn="Ramp", mod_names=mod_names)
axis(1, at=c(1,2), labels=mod_names, cex.axis=2)
mtext("Relative Error", side=2, line=2, cex=1.3, outer=TRUE)

linf_sites <- y_sites <- matrix(NA, nrow=n_i, ncol=64)
for(i in 1:64){
	sim <- readRDS(file.path(modpath_space, "F_Constant", i, "spatial_sim.rds"))
	linf_sites[,i] <- sim$linf_i
	y_sites[,i] <- sim$y_i
}
par(mfrow=c(8,8), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:64){
	plot(y_sites[,i], linf_sites[,i],axes=F, ann=F, ylim=c(min(linf_sites), max(linf_sites)), pch=19)
	box(col="gray")
	abline(h=linf, col="red", lty=2)
}

setwd(data_dir)
source("R_functions\\functions.R")

plotFIT(compare_quant=c( "ML", "R", "F", "D"), compare_type="base_values", scenario_name="No spatial structure", modpath=file.path(modpath_nospace, "F_Constant"), iter=1)
dev.new()
plotFIT(compare_quant=c( "ML", "R", "F", "D"), compare_type="base_values", scenario_name="Spatial structure", modpath=file.path(modpath_space, "F_Constant"), iter=1)


############ Part 2 ###############
### QUESTION
### Is it more important to account for temporal variation in fishing mortality, temporal variation in recruitment, or spatial variation in growth for estimates of stock status?

### PERFORMANCE MEASURE
#### the ability to estimate the mean SPR across sites with spatial variation

### METHOD STRUCTURE
## OM1 = spatial differences in growth, constant recruitment with process error, constant/endogenous fishing mortality with process error
## EM1 = estimate annual F only (turn off random effects on R and linf)
## EM2 = estimate variable R only (estimate only 1 F, turn off random effects on linf)
## EM3 = estimate variable linf only (estimate only 1 F, turn off random effects on R)
## EM4 = estimate annual F and variable R (turn off random effects on linf)
## EM5 = estimate annual F and variable linf (turn off random effects on R)
## EM6 = estimate variable R and linf (estimate only 1 F)
rm(list=ls())

# data_dir <- "C:\\Git_Projects\\2016_Spatio-temporal_models\\Project"
data_dir <- "F:\\Merrill\\Git_Projects\\2016_Spatio-temporal_models\\Project"
setwd(data_dir)
source("R_functions\\functions.R")

require(geoR)
require(RandomFields)
require(ape)
require(RColorBrewer)
require(spdep)
require(reshape2)
require(TMB)

## locate and compile executable
run_exe <- file.path(data_dir, "inst", "executables")

Version <- "lb_statespace_v6"
setwd(run_exe)
dyn.unload( paste0(Version,".dll") )
file.remove( paste(Version,c(".dll",".o"),sep="") )
compile(paste0(Version, ".cpp"))

## directories for results
estimate_variability <- file.path(data_dir, "estimate_variability")
dir.create(estimate_variability, showWarnings=FALSE) 

unlink(estimate_variability, TRUE)
dir.create(estimate_variability, showWarnings=FALSE) 

timeF <- file.path(estimate_variability, "timeF")
dir.create(timeF, showWarnings=FALSE)

timeR <- file.path(estimate_variability, "timeR")
dir.create(timeR, showWarnings=FALSE)

spaceL <- file.path(estimate_variability, "spaceL")
dir.create(spaceL, showWarnings=FALSE)

timeF_timeR <- file.path(estimate_variability, "timeF_timeR")
dir.create(timeF_timeR, showWarnings=FALSE)

timeF_spaceL <- file.path(estimate_variability, "timeF_spaceL")
dir.create(timeF_spaceL, showWarnings=FALSE)

timeR_spaceL <- file.path(estimate_variability, "timeR_spaceL")
dir.create(timeR_spaceL, showWarnings=FALSE)

## settings - general
dat_scenario <- "Poor_Comp_LC"
Nyears <- 20
comp_sample <- 100 
SigmaR <- 0.6
seed <- 123
Nyears_comp <- 10
simdir <- "pool_space"
sens_param <-FALSE
sens_val <- FALSE
OMselex <- "asymptotic"
EMselex <- "asymptotic"
nburn <- 20
est_sigma <- FALSE
rerun=FALSE
runLBSPR=FALSE
REML=FALSE

## settings - space
### simulation of 1D spatial process on growth - from HW5
n_i <- 15 ##  number of sites
Scale <- 2
Sigma2 <- 0.1
SD_spatial <- 0.1
linf <- 36.2 
beta_y <- 0.05

setwd(data_dir)
source("R_functions\\functions.R")

## model path
modpath_timeF <- file.path(estimate_variability, "timeF")
modpath_timeR <- file.path(estimate_variability, "timeR")
modpath_spaceL <- file.path(estimate_variability, "spaceL")
modpath_timeF_timeR <- file.path(estimate_variability, "timeF_timeR")
modpath_timeF_spaceL <- file.path(estimate_variability, "timeF_spaceL")
modpath_timeR_spaceL <- file.path(estimate_variability, "timeR_spaceL")

## life history - truth without spatial structure - deterministic across iterations
lh_nospace <- create_lh_list(lh="Siganus_sutor", selex="asymptotic")
## data inputs - same across scenarios, assumes no spatial differences in growth
dat_input <- create_inputs(param=FALSE, val=FALSE, lh_dat=lh_nospace)

niter <- 1

## generate data - spatial structure of linf in true population, Linf pooled (not accounting for spatial growth in EM)
generateData(modpath=modpath_timeF, itervec=1:niter, spatial=TRUE, Fdynamics="Constant", LType=0)
	check <- readRDS(file.path(modpath_timeF, "F_constant", "1", "True.rds"))
## run model - time-varying F
runModel(modpath=modpath_timeF, itervec=1:niter, data_input=dat_input, Fdynamics="Constant", RecType=1, FType=0, LType=1, site=1)

# 	setwd(data_dir)
# 	source("R_functions\\functions.R")


# 	

# 	## compare SPR at each site with Linf at that site
# 	par(mfrow=c(2,1))
# 	plot(x=spatial_sim$y_i, y=spatial_sim$linf_i)
# 	plot(x=1:length(DataList_site), y=sapply(1:length(DataList_site), function(x) DataList_site[[x]]$SPR), ylim=c(0,1))





### explore real data
SSraw <- read.csv(file.path(data_dir, "Siganus_sutor_data.csv"), header=TRUE)

## Siganus sutor by region
SSnew <- SSraw
SSnew$Landing <- tolower(SSraw$Landing)
SS_locs <- unique(SSnew$Landing)
years <- min(SSnew$Year):max(SSnew$Year)
SSdat <- createDataList(sp="Siganus_sutor")

SS2014 <- SSnew[which(SSnew$Year==2014),]
SS_locs2014 <- unique(SS2014$Landing)

lcomp_site <- list()
for(i in 1:length(SS_locs2014)){
	sub_site <- SS2014[which(SS2014$Landing==SS_locs2014[i]),]
	lengths <- as.numeric(unlist(sub_site[,which(colnames(sub_site)=="X1"):which(colnames(sub_site)=="X32")]))
	lcomp_site[[SS_locs2014[i]]] <- lengths[-which(is.na(lengths))]
}

par(mfrow=c(5,3))
breaks <- seq(0,40,by=2)
lengths_all <- as.numeric(unlist(SS2014[,which(colnames(SS2014)=="X1"):which(colnames(SS2014)=="X32")]))
x <- hist(lengths_all[-which(is.na(lengths_all))], plot=FALSE, breaks=breaks)
plot(x=1, y=1, type="n", xlim=c(min(breaks),max(breaks)), ylim=c(0,0.5), xaxs="i", yaxs="i", xaxt="n", yaxt="n")
polygon(x=c(breaks,rev(breaks)), y=c(0,x$counts/sum(x$counts),rep(0,length(x$density)+1)), col="blue")

for(i in 1:length(SS_locs2014)){
	x <- hist(lcomp_site[[i]], plot=FALSE, breaks=breaks)
plot(x=1, y=1, type="n", xlim=c(min(breaks),max(breaks)), ylim=c(0,0.5), xaxs="i", yaxs="i", xaxt="n", yaxt="n")
polygon(x=c(breaks,rev(breaks)), y=c(0,x$counts/sum(x$counts),rep(0,length(x$density)+1)), col="blue")

}

## mean length by year and site
ml_yr_site <- matrix(NA, nrow=length(SS_locs), ncol=length(years))
colnames(ml_yr_site) <- years
rownames(ml_yr_site) <- SS_locs
for(i in 1:length(SS_locs)){
	sub_site <- SSnew[which(SSnew$Landing==SS_locs[i]),]
	yrs_site <- unique(sub_site$Year)
	for(j in 1:length(yrs_site)){
		sub_yr <- sub_site[which(sub_site$Year==yrs_site[j]),]
		lengths <- as.numeric(unlist(sub_yr[,which(colnames(sub_yr)=="X1"):which(colnames(sub_yr)=="X32")]))
		ml_yr_site[which(rownames(ml_yr_site)==SS_locs[i]), which(colnames(ml_yr_site)==yrs_site[j])] <- mean(lengths, na.rm=TRUE)
	}
}

## mean length by year
ml_yr <- rep(NA, length(years))
names(ml_yr) <- years
for(j in 1:length(years)){
	sub_yr <- SSnew[which(SSnew$Year==years[j]),]
	lengths <- as.numeric(unlist(sub_yr[,which(colnames(sub_yr)=="X1"):which(colnames(sub_yr)=="X32")]))
	ml_yr[which(names(ml_yr)==years[j])] <- mean(lengths, na.rm=TRUE)
}

png("Mean_length_by_site_year.png", res=200, units="in", width=12, height=9)
par(mfrow=c(4,6), mar=c(0,0,0,0), omi=c(1,1,1,1))
for(i in 1:length(SS_locs)){
	plot(x=years, y=ml_yr_site[i,], pch=19, cex=2, xaxs="i", yaxs="i", ylim=c(0, max(ml_yr_site, na.rm=TRUE)*1.2), xaxt="n", yaxt="n", xpd=NA, xlab="", ylab="")
	if(i %in% 19:24) axis(1, cex.axis=1.3)
	if(i %in% c(1,7,13,19)) axis(2, las=2, cex.axis=1.3)
	mtext(side=3, rownames(ml_yr_site)[i], font=2, cex=1.3, line=-2)
	lines(x=names(SSdat$meanlen), y=SSdat$meanlen, col="red", pch=17, type="o")
	if(i==length(SS_locs)) legend("bottomright", legend=c("site-specific", "pooled"), pch=c(19,17), cex=1.2, col=c("black", "red"))
}
mtext("Year", outer=TRUE, side=1, line=3)
mtext("Mean length in catch (cm)", side=2, line=3, outer=TRUE)
dev.off()
