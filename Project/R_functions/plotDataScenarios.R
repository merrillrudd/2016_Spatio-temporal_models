plotDataScenarios <- function(){

	maxyrs <- 20
xlim <- c(1, maxyrs)
ylim <- c(0, 15)
    yrs <- 1:maxyrs
    ryrs <- rev(yrs)
    index <- rep(c(1,1,0,0), maxyrs/4)
    alt_yrs <- rev(ryrs[which(index==1)])

par(mfrow=c(1,2), mar=c(4,2,2,7), omi=c(0,0,0,1))
plot(x=xlim[1]:xlim[2], y=seq(ylim[1], ylim[2], length=length(xlim[1]:xlim[2])), 
  type="n", ann=F, xaxs="i", yaxs="i", axes=F)
text(x=10, y=14, "Rich", font=2, cex=1.2)
lines(x=yrs, y=rep(13, maxyrs), col="blue", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=13, "Catch", xpd=NA, cex=1.2)
lines(x=yrs, y=rep(12, maxyrs), col="red", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=12, "Index", xpd=NA, cex=1.2)
points(x=yrs, y=rep(11, maxyrs), col="goldenrod", pch=19, cex=1000/200, xpd=NA)
text(x=(maxyrs+4), y=11, "Length Comp", xpd=NA, cex=1.2)
text(x=(maxyrs+4.5), y=10.4, "ESS = 500", xpd=NA, cex=1.2)

text(x=10, y=9, "Moderate (Full)", font=2, cex=1.2)
lines(x=yrs, y=rep(8, maxyrs), col="blue", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=8, "Catch", xpd=NA, cex=1.2)
lines(x=yrs, y=rep(7, maxyrs), col="red", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=7, "Index", xpd=NA, cex=1.2)
points(x=yrs, y=rep(6, maxyrs), col="goldenrod", pch=19, cex=3, xpd=NA)
text(x=(maxyrs+4), y=6, "Length Comp", xpd=NA, cex=1.2)
text(x=(maxyrs+4), y=5.4, "ESS = 50", xpd=NA, cex=1.2)

text(x=10, y=4, "Moderate (Sample)", font=2, cex=1.2)
points(x=alt_yrs, y=rep(3, length(alt_yrs)), col="blue", xpd=NA, pch=19, lwd=10)
text(x=(maxyrs+2.2), y=3, "Catch", xpd=NA, cex=1.2)
points(x=alt_yrs, y=rep(2, length(alt_yrs)), col="red", xpd=NA, pch=19, lwd=10)
text(x=(maxyrs+2.2), y=2, "Index", xpd=NA, cex=1.2)
points(x=alt_yrs, y=rep(1, length(alt_yrs)), col="goldenrod", pch=19, cex=3, xpd=NA)
text(x=(maxyrs+4), y=1, "Length Comp", xpd=NA, cex=1.2)
text(x=(maxyrs+4), y=0.4, "ESS = 50", xpd=NA, cex=1.2)
axis(1, at=c(1,seq(5,20,by=5)))
mtext("Year", side=1, line=2, cex=1.2)


# par(mfrow=c(1,1), omi=c(0,0,0,2))
plot(x=xlim[1]:xlim[2], y=seq(ylim[1], ylim[2], length=length(xlim[1]:xlim[2])), 
  type="n", ann=F, xaxs="i", yaxs="i", axes=F)
text(x=10, y=14, "Poor (Index + Comp)", font=2, cex=1.2)
# lines(x=yrs, y=rep(13, maxyrs), col="blue", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=13, "Catch", xpd=NA, cex=1.2)
lines(x=11:20, y=rep(12, 10), col="red", lwd=10, xpd=NA)
text(x=(maxyrs+2.2), y=12, "Index", xpd=NA, cex=1.2)
points(x=maxyrs, y=11, col="goldenrod", pch=19, cex=3/2, xpd=NA)
text(x=(maxyrs+4), y=11, "Length Comp", xpd=NA, cex=1.2)
text(x=(maxyrs+4.5), y=10.4, "ESS = 10", xpd=NA, cex=1.2)

# text(x=10, y=9, "Poor (Comp + Relative)", font=2, cex=1.2)
# points(x=c(maxyrs-19,maxyrs), y=rep(8, 2), col="blue", xpd=NA, pch=19, lwd=5)
# text(x=(maxyrs+2.2), y=8, "Catch", xpd=NA, cex=1.2)
# points(x=c(maxyrs-19,maxyrs), y=rep(7, 2), col="red", xpd=NA, pch=19, lwd=5)
# text(x=(maxyrs+2.2), y=7, "Index", xpd=NA, cex=1.2)
# points(x=maxyrs, y=6, col="goldenrod", pch=19, cex=3/2, xpd=NA)
# text(x=(maxyrs+4), y=6, "Length Comp", xpd=NA, cex=1.2)
# text(x=(maxyrs+4), y=5.4, "ESS = 10", xpd=NA, cex=1.2)

text(x=10, y=4, "Poor (Comp only)", font=2, cex=1.2)
text(x=(maxyrs+2.2), y=3, "Catch", xpd=NA, cex=1.2)
text(x=(maxyrs+2.2), y=2, "Index", xpd=NA, cex=1.2)
points(x=1:20, y=rep(1, 20), col=c(rep("gray",19),"goldenrod"), pch=19, cex=3/2, xpd=NA)
text(x=(maxyrs+4), y=1, "Length Comp", xpd=NA, cex=1.2)
text(x=(maxyrs+4), y=0.4, "ESS = 10", xpd=NA, cex=1.2)

axis(1, at=c(1,seq(5,20,by=5)))
mtext("Year", side=1, line=2, cex=1.2)


}