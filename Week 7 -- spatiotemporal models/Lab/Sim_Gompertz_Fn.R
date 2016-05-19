Sim_Gompertz_Fn <-
function( n_years, n_stations=100, phi=NULL, SpatialScale=0.1, SD_O=0.5, SD_E=0.2, SD_extra=0.1, rho=0.8, logMeanDens=1, Loc=NULL ){
  # Parameters
  if( is.null(phi) ) phi = rnorm(1, mean=0, sd=1)
  alpha = logMeanDens * (1-rho)

  # Spatial model -- implicitly assuming rho between -1 and 1
## Loc - simplest random sampling
## RMgauss from RandomFields package - creates an object with a variance and scale (fed in), then use RFsimulate (2D random field) to simulate value of omega at locations that we have randomly generated from a process with a certain standard deviation and scale
  if( is.null(Loc) ) Loc = cbind( "x"=runif(n_stations, min=0,max=1), "y"=runif(n_stations, min=0,max=1) )
  model_O <- RMgauss(var=SD_O^2, scale=SpatialScale)
  model_E <- RMgauss(var=SD_E^2, scale=SpatialScale)

  # Simulate Omega
  Omega = RFsimulate(model = model_O, x=Loc[,'x'], y=Loc[,'y'])@data[,1]

  # Simulate Epsilon - sweep upstream to downstream by year to generate epsilon
  Epsilon = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    Epsilon[,t] = RFsimulate(model=model_E, x=Loc[,'x'], y=Loc[,'y'])@data[,1]
  }

  # Calculate Psi - same size as epsilon --- equilibrium in year 1 (need omega in spatial model)
  Theta = array(NA, dim=c(n_stations,n_years))
  for(t in 1:n_years){
    if(t==1) Theta[,t] = phi + Epsilon[,t] + (alpha + Omega)/(1-rho)
    if(t>=2) Theta[,t] = rho * Theta[,t-1] + alpha + Omega + Epsilon[,t]
  }

  # Simulate data
  DF = NULL
  for(s in 1:n_stations){
  for(t in 1:n_years){
    Tmp = c("Site"=s, "Year"=t, "Simulated_example"=rpois(1,lambda=exp(Theta[s,t]+SD_extra*rnorm(1)))) ## includes overdispersion (SD_extra -- same as doing rnorm(1, 0, SD_extra))
    DF = rbind(DF, Tmp)
  }}
  ## record 'pretend' latitude and longitude of each site
  DF = cbind( DF, 'Longitude'=Loc[DF[,'Site'],1], 'Latitude'=Loc[DF[,'Site'],2] )
  DF = data.frame(DF, row.names=NULL)

  # Return stuff
  Sim_List = list("DF"=DF, "phi"=phi, "Loc"=Loc, "Omega"=Omega, "Epsilon"=Epsilon, "Theta"=Theta)
  return(Sim_List)
}
