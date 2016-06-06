## spatial process in Linf - incorporates anabolism and catabolism
spatialgrowth_sim <- function(n_i, Scale, Sigma2, SD_spatial, linf, beta_y){
	require(RandomFields)

	## sample locations
	lat_min <- -4
	lat_max <- -1
	y_i <- runif(n=n_i, min=lat_min, max=lat_max)

	## simulate spatial process
	RMmodel <- RMgauss(var=SD_spatial^2, scale=Scale)
	linf_i <- linf * exp(RFsimulate(model=RMmodel, x=rep(0, n_i), y=y_i)@data[,1] - Sigma2/2) * exp( beta_y*(y_i - mean(c(lat_min, lat_max))))

	## resimulate to ensure Linf > length at 50% maturity from literature
	while(any(linf_i < 20.2)){
		RMmodel <- RMgauss(var=SD_spatial^2, scale=Scale)
		linf_i <- linf * exp(RFsimulate(model=RMmodel, x=rep(0, n_i), y=y_i)@data[,1] - Sigma2/2) * exp( beta_y*(y_i - mean(y_i)))
	}
	# plot( y=linf_i, x=y_i )


	df <- data.frame(linf_i=linf_i, y_i=y_i)
	return(df)
}