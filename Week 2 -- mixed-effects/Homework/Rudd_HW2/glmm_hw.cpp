#include <TMB.hpp>


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( model ); // choose estimation model
  DATA_INTEGER( tobs ); // total number of observations
  DATA_INTEGER( nsites ); // number of sites
  DATA_VECTOR( y_j ); // observations
  DATA_IVECTOR( site_j ); // sites

  // Parameters
  PARAMETER( mu ); // log-mean for expected counts
  PARAMETER( logsig_s );  // standard deivation of among-site variability
  PARAMETER( logsig_y ); // standard deivation of overdispersion

  // Random effects
  PARAMETER_VECTOR( tau ); // density at each site
  PARAMETER_VECTOR( epsilon ); // observation at each site

  // transformations
  Type sig_s = exp(logsig_s);
  Type sig_y = exp(logsig_y);

  // joint likelihood
  Type jnll = 0;
  vector<Type> jnll_comp(3);

  // probability of random effects
  // density at each site
  jnll_comp(0) = 0;
  if(model==2 | model==4){
    for(int i=0; i<nsites; i++){
      jnll_comp(0) -= dnorm(tau(i), Type(0.0), sig_s, true);
    }
  }
  // overdispersion
  jnll_comp(1) = 0; 
  if(model==3 | model==4){
    for(int j=0; j<tobs; j++){
      jnll_comp(1) -= dnorm(epsilon(j), Type(0.0), sig_y, true);
    }
  }

  // prediction
  vector<Type> pred_j( tobs );
  for(int j=0; j<tobs; j++){
    if(model==1) pred_j(j) = exp(mu);
    if(model==2) pred_j(j) = exp(mu + tau(site_j(j)));
    if(model==3) pred_j(j) = exp(mu + epsilon(j));
    if(model==4) pred_j(j) = exp(mu + tau(site_j(j)) + epsilon(j));

    // if(model==2) pred_j(j) = exp(mu + tau(site_j(j)) - pow(sig_s,2)/Type(2));
    // if(model==3) pred_j(j) = exp(mu + epsilon(j) - pow(sig_y,2)/Type(2));
    // if(model==4) pred_j(j) = exp(mu + tau(site_j(j)) - pow(sig_s,2)/Type(2) + epsilon(j) - pow(sig_y,2)/Type(2));
  }

  jnll_comp(2) = 0;
  for(int j=0; j<tobs; j++){
    jnll_comp(2) -= dpois(y_j(j), pred_j(j), true);
  }

  jnll = sum(jnll_comp);

  REPORT(jnll);
  REPORT(jnll_comp);
  REPORT(model);

  REPORT(mu);
  REPORT(sig_s);
  REPORT(sig_y);
  REPORT(tau);
  REPORT(epsilon);

  REPORT(pred_j);
  ADREPORT(pred_j);

  return(jnll);

}
