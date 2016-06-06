createDataList <- function(sp){
    data <- read.csv(file.path(data_dir, paste0(sp, "_data.csv")), header=TRUE)
    meltDat <- melt(data, id.vars=c("Date", "Year", "Landing", "Management.type", "Group", 
    	"Fishgear","X.Fishers", "Catch.category", "Family", "Name", "Catch.category",
    	"Commercial.noncomm."))
    meltDat$value <- as.numeric(meltDat$value)
    dat_new <- meltDat[-which(is.na(meltDat$value)),-which(colnames(meltDat)=="variable")]    


    years <- unique(dat_new$Year)[order(unique(dat_new$Year))]
    lbins <- 0:(max(dat_new$value)*1.5)    

    lcomp1 <- matrix(0, nrow=length(years), ncol=length(lbins))
    colnames(lcomp1) <- lbins
    rownames(lcomp1) <- years   
    lcomp_p <- lcomp1 

    catch <- cpue <- vector(length=length(years))
    cpue_day_out <- list()    

    for(y in 1:length(years)){
    	sub <- dat_new[which(dat_new$Year==years[y]),]    

    	## annual catch
    	catch[y] <- nrow(sub)    

    	## annual cpue - average of daily/group cpue
    	cpue_sub <- rep(NA, length(unique(sub$Date)))
        ngroups <- rep(NA, length(unique(sub$Date)))
    	for(d in 1:length(unique(sub$Date))){
    		subd <- sub[which(sub$Date==unique(sub$Date)[d]),]    
            ngroups <- length(unique(subd$Group)[which(is.na(unique(subd$Group))==FALSE)])
    	    if(ngroups==0) cpue_sub[d] <- NA
            if(ngroups!=0) cpue_sub[d] <- nrow(subd)/ngroups
            # rm(ngroups)
            # rm(subd)
        }
    	if(all(is.na(cpue_sub))) cpue[y] <- NA
    	if(all(is.na(cpue_sub)==FALSE)) cpue[y] <- mean(cpue_sub, na.rm=TRUE)
        cpue_day_out[[y]] <- cpue_sub
    	rm(cpue_sub)    

        for(l in 1:length(lbins)){
        	lcomp1[y,l] <- length(which(floor(sub$value)==lbins[l]))
        }
    }    

    names(cpue) <- names(catch) <- years
    cpue <- cpue[-c(which(is.na(cpue)),which(cpue==0))]  

    meanlen <- rep(NA, nrow(lcomp1))
    names(meanlen) <- rownames(lcomp1)
    for(i in 1:nrow(lcomp1)){
        lcomp_p[i,] <- lcomp1[i,]/sum(lcomp1[i,])
        meanlen[i] <- sum(sapply(1:ncol(lcomp1), function(x) lcomp1[i,x]*as.numeric(colnames(lcomp1)[x])))/sum(lcomp1[i,])
    }  

    DataList <- NULL
    DataList$I_t <- cpue
    DataList$C_t <- catch
    DataList$LF <- lcomp1
    DataList$LFprop <- lcomp_p
    DataList$years <- years
    DataList$lbins <- lbins
    DataList$meanlen <- meanlen
    return(DataList)

}