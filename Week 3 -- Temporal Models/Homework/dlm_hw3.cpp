
// State-space dynamic linear model
#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
// data:
DATA_INTEGER(model);
DATA_VECTOR(y);

// parameters:
PARAMETER(a); // population growth rate parameter
PARAMETER(b);
PARAMETER(log_sigma_proc); // log(process SD)
PARAMETER(log_sigma_obs); // log(observation SD)
PARAMETER_VECTOR(u); // unobserved state


// transform the parameters
Type sigma_proc = exp(log_sigma_proc);
Type sigma_obs = exp(log_sigma_obs);

// reports
ADREPORT(sigma_proc)
ADREPORT(sigma_obs)

int n = y.size(); // get time series length

Type nll = 0.0; // initialize negative log likelihood

// process model:
vector<Type> expected(n);
expected(0) = u(0);
for(int i = 1; i < n; i++){
	//linear model
  if(model==0) expected(i) = a*u[i-1] ;    // calculated mean of the process, 'expected' keeping track of predicted value, state at time t-1, recast every time the loop goes through 
  if(model==1) expected(i) = a + b*u[i-1];
  nll -= dnorm(u[i], expected(i), sigma_proc, true); //likelihood for random effects
}

// observation model:
for(int i = 0; i < n; i++){
  nll -= dnorm(y[i], u[i], sigma_obs, true); //likelihood for observations
}

return nll;
}
