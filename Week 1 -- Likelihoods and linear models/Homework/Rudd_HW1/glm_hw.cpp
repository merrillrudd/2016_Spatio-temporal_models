#include <TMB.hpp>

// dlnorm
template<class Type>
Type dlognorm(Type x, Type meanlog, Type sdlog, int give_log=0){
  Type logres = dnorm( log(x), meanlog, sdlog, true) - log(x);
  if(give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data
  DATA_INTEGER( like ); // define likelihood type, 1==delta lognormal, 2==delta gamma
  DATA_VECTOR( y_i ); // observations
  DATA_MATRIX( X_ij ); // covariate design matrix
  DATA_VECTOR( include ); //0== include in NLL, 1== exclude from NLL

  // Parameters
  PARAMETER_VECTOR( b_j ); // betas to generate expected values
  PARAMETER_VECTOR( theta_z ); // variances

  // Transformations
  Type zero_prob = 1 / (1 + exp(-theta_z(0)));
  Type sd = exp(theta_z(1)); //standard deviation (lognormal), scale parameter theta (gamma)
  int n_data = y_i.size();

  Type jnll = 0;
  Type pred_jnll = 0;
  vector<Type> jnll_i(n_data);


  // linear predictor
  vector<Type> logpred_i( n_data );
  logpred_i = X_ij * b_j; 


  // Delta lognormal
  if(like==1){
  	for( int i=0; i<n_data; i++){
      if(y_i(i)==0) jnll_i(i) -= log( zero_prob );
      if(y_i(i)!=0) jnll_i(i) -= log( 1-zero_prob ) + dlognorm( y_i(i), logpred_i(i), sd, true );
    
      // Running counter
      if( include(i)==0 ) jnll += jnll_i(i);
      if( include(i)==1 ) pred_jnll += jnll_i(i);
    }
  }

  // Delta gamma
  if(like==2){
    for(int i=0; i<n_data; i++){
    	if(y_i(i)==0) jnll_i(i) -= log( zero_prob );
    	if(y_i(i)!=0) jnll_i(i) -= log( 1-zero_prob ) + dgamma( y_i(i), pow(sd,-2), exp(logpred_i(i))*pow(sd,2), true );

        // Running counter
        if( include(i)==0 ) jnll += jnll_i(i);
        if( include(i)==1 ) pred_jnll += jnll_i(i);
    }
  }



  
  // Reporting
  REPORT( zero_prob );
  REPORT( sd );
  REPORT( logpred_i );
  REPORT( b_j );
  REPORT( pred_jnll );
  REPORT( jnll_i );
  return jnll;
}
