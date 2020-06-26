// An attempt to translate a jags model into stan. I successfully ran stan
// models for the leatherback juvenile survival analysis so I should be able
// to do this also. Started 2020-06-26



data {
  
  int<lower=1> n_months;
  int<lower=1> n_years;  // n.steps in the original jags code
  //int<lower=1> n_timeseries;
  int<lower=1> period_1;

  int<lower=1> m[12 * n_years];
  
  int<lower = 0> N_obs;  // # observed data points
  int<lower = 0> N_miss;  // # missing data points

  int<lower = 1, upper = N_obs + N_miss> idx_obs[N_obs];  // index of observed data
  int<lower = 1, upper = N_obs + N_miss> idx_miss[N_miss]; // index of missed data
  
  real<lower = 0.0> y_obs[N_obs];
  
  real c_const[period_1];
  
}

transformed data {
  int<lower = 0> N = N_obs + N_miss; 
  
}

parameters {
  real<lower = 0.0> y_miss[N_miss];
  
  real<lower = 0> predX0;
  real predX[N];
  real<lower = 0.0> X[N];
  
  real<lower=0> c[n_months];
  real beta_cos;
  real beta_sin;
  real<lower=0> sigma_X;
  real<lower=0> sigma_y;
  
}

transformed parameters {
  
  real<lower = 0.0> y[N];
  y[idx_obs] = y_obs;
  y[idx_miss] = y_miss;

}	

model {
  
  // Initial state space
  predX0 ~ normal(8, 10);
  X[1] ~ normal(beta_cos * cos(c_const[m[1]]) + beta_sin * sin(c_const[m[1]]) + predX0, sigma_X);
  
  // observation
  y[1] ~ normal(X[1], sigma_y);
                
  for (t in 2:N){
    // State space
    X[t] ~ normal(beta_cos * cos(c_const[m[t]]) + beta_sin * sin(c_const[m[t]]) + X[t-1], sigma_X);
    
    // observation
    y[t] ~ normal(X[t], sigma_y);
    
  }
       
  sigma_y ~ gamma(2, 0.5);
  beta_cos ~ normal(0, 1);
  beta_sin ~ normal(0, 1);

  // sigma.X is the SD of the process (X)
  sigma_X ~ gamma(2, 0.5);
  
}

generated quantities {
  // log_lik is for use with the loo package
  vector[N_obs] log_lik;
  for(i in 1:N_obs) {
  	log_lik[i] = lognormal_lpdf(y[i] | X[i], sigma_y);
  }
}

