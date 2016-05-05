
#include <TMB.hpp>

template<class Type>
Type dmvnorm( vector<Type> x, matrix<Type> Q, int give_log=0 ){
  int n_x = x.size();
  Type logres = 0;
  vector<Type> Tmp_x(n_x);
  Type logdet_Q = logdet( Q );
  Tmp_x =  Q * x.matrix();
  logres += ( Type(0.5)*logdet_Q );
  logres += Type(-0.5) * (x * Tmp_x).sum();  //
  if (give_log) return logres; else return exp(logres);
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_INTEGER(Option);

  // Data
  DATA_VECTOR( a_i );
  DATA_VECTOR( l_i );
  DATA_VECTOR( loc_i );

  // Parameters
  PARAMETER( beta0 );
  PARAMETER( betaS );
  PARAMETER( l0 );
  PARAMETER( vbk );
  PARAMETER( ln_sigma2 );
  PARAMETER( logit_rho );

  // Random effects
  PARAMETER_VECTOR( epsilon_i );

  // Objective function and transformations
  int n_i = a_i.size();   // number of samples
  vector<Type> jnll_comp(2);
  jnll_comp.setZero();
  Type sigma2 = exp(ln_sigma2);
  Type rho = 1 / (1 + exp(-logit_rho));

  // Probability of random effects - linf
  using namespace density;

  if(Option==1){
    vector<Type> dist_i(n_i);
    jnll_comp(1) -= dnorm( epsilon_i(0), Type(0.0), pow(sigma2,0.5), true );
    for(int i=1; i<n_i; i++){
      dist_i(i) = loc_i(i)-loc_i(i-1);
      jnll_comp(1) -= dnorm( epsilon_i(i), pow(rho,dist_i(i))*epsilon_i(i-1), pow(sigma2*(1-pow(rho,2*dist_i(i))),0.5), true );
    }
  }

  if(Option==2){
    // TMB multi-variate normal function
    matrix<Type> Cov_ss(n_i,n_i);
    Type x;
    for(int s1=0; s1<n_i; s1++){
    for(int s2=s1; s2<n_i; s2++){
      Cov_ss(s1,s2) = sigma2 * pow( rho, abs(loc_i(s2)-loc_i(s1)) );
      if(s1!=s2) Cov_ss(s2,s1) = Cov_ss(s1,s2);
    }}
    jnll_comp(1) += MVNORM(Cov_ss)( epsilon_i );
      REPORT( Cov_ss );

  }

  // predicted linf by observation
  vector<Type> linf_i(n_i);
  for(int i=0; i<n_i; i++){
    linf_i(i) = exp(beta0 + betaS*loc_i(i) + epsilon_i(i));
  }

  // predicted length at age - von Bertalanffy
  vector<Type> pred_l_i(n_i);
  for(int i=0; i<n_i; i++){
    pred_l_i(i) = linf_i(i) - (linf_i(i) - l0) * exp(-vbk*a_i(i));
  }

  // probability of data conditional on random effects
  for(int i=0; i<n_i; i++){
    jnll_comp(0) -= dnorm(log(l_i(i)), log(pred_l_i(i)), pow(sigma2,0.5), true);
  }

  // Reporting
  Type jnll = jnll_comp.sum();
  REPORT( jnll_comp );
  REPORT( jnll );
  REPORT( sigma2 );
  REPORT( rho );
  REPORT( linf_i );
  REPORT( pred_l_i );

  return jnll;
}

