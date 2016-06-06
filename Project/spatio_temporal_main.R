rm(list=ls())

library(reshape2)
library(TMB)

data_dir <- work_dir <- "C:\\Git_Projects\\kenyan_reef_fish"

setwd(work_dir)
source("R_functions\\functions.R")


###########################
## Data
###########################

  sp_vec <- c("Siganus_sutor", "Lethrinus_lentjan", "Leptoscarus_vaigiensis")
  Lmat <- c("Siganus_sutor"=20.2, "Lethrinus_lentjan"=20.3, "Leptoscarus_vaigiensis"=15.1)

  SSdat <- createDataList(sp=sp_vec[1])
  LLdat <- createDataList(sp=sp_vec[2])
  LVdat <- createDataList(sp=sp_vec[3])

  SSraw <- read.csv(file.path(data_dir, paste0(sp_vec[1], "_data.csv")), header=TRUE)
  LLraw <- read.csv(file.path(data_dir, paste0(sp_vec[2], "_data.csv")), header=TRUE)
  LVraw <- read.csv(file.path(data_dir, paste0(sp_vec[3], "_data.csv")), header=TRUE)


## mean length of 3 species
par(mfrow=c(1,1), omi=c(1,2,1,1))
 plot(x=names(SSdat$meanlen), y=SSdat$meanlen/Lmat[1], ylim=c(0, max(SSdat$meanlen)/Lmat[1]*1.2), type="o", pch=19, cex=1.2, xaxs="i", yaxs="i", xlim=c(min(as.numeric(names(SSdat$meanlen))), max(as.numeric(names(SSdat$meanlen)))+1), xlab="", ylab="", lwd=2)
 lines(x=names(LLdat$meanlen), y=LLdat$meanlen/Lmat[2], type="o", pch=19, cex=1.2, col="steelblue", lwd=2)
 lines(x=names(LVdat$meanlen), y=LVdat$meanlen/Lmat[3], type="o", pch=19, cex=1.2, col="tomato", lwd=2)
 mtext(side=1, "Year", cex=1.5, line=3)
 mtext(side=2, "Mean length in catch\nrelative to length at maturity", cex=1.5, line=3)
 legend("bottomright", legend=c("Siganus sutor", "Lethrinus lentjan", "Leptoscarus vaigiensis"), col=c("black", "steelblue", "tomato"), lwd=2, lty=1, pch=19, cex=1.5)
 abline(h=1, lty=2)

## Siganus sutor by region
SSnew <- SSraw
SSnew$Landing <- tolower(SSraw$Landing)
SS_locs <- unique(SSnew$Landing)

SS_reg <- c("shimoni"="south", "mwaepe"="south", "chale"="south", "mvuleni"="south", 
	"nyali"="south", "msanakani"="south", "gazi"="south", "kuruwitu"="south", "kijangwani"="south", "kanamai"="south")

### plot mean length over time by landings site
### map sites -- some inland? which are southern coast? which north?
