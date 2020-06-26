data {
	int<lower=1> Y;          // The total number of years in the time series that can be used 
	real<lower=0> N[Y+3];   // # of neophytes 9 years ago. although integers, modeled as a real
	real<lower=0> H[Y];    // # of observed hatchlings. although integers, modeled as a real
	
}

parameters {
	real<lower=0, upper = 10> sigma;       // standard deviation of log-normal hatchling counts
	real<lower=0, upper = 1.0> phi[Y];        // annual average survival rates for juveniles	
	simplex[4] q[Y+3];
}

transformed parameters {
	real<lower=0> mu[Y];     // the log-mean of the # hatchlings 
	real<lower=0> alpha = 1.0;   // dirichlet parameters

	for (i in 1:Y){
		mu[i] = log(q[i, 1] * N[i] * pow(phi[i],-9) + 
		             q[i+1, 2] * N[i+1] * pow(phi[i],-10) +
		             q[i+2, 3] * N[i+2] * pow(phi[i],-11) + 
		             q[i+3, 4] * N[i+3] * pow(phi[i],-12)); 

	}	

	
}

model
{
	// Priors     	
	sigma ~ uniform(0 , 10);
	//phi ~ beta(1,1);

	for (i in 1:Y){
		phi[i] ~ beta(1,1);
		H[i] ~ lognormal(mu[i], sigma);
		q[i] ~ dirichlet(rep_vector(alpha, 4));
		

	}
  
}

generated quantities {
  // log_lik is for use with the loo package
  vector[Y] log_lik;
  for(i in 1:Y) {
  	log_lik[i] = lognormal_lpdf(H[i] | mu[i], sigma);
  }
}

