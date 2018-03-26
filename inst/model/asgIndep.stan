data {
  int<lower=1> n; //total number of obs over all locations and seasons
  int<lower=1> num[n]; //count of total patients
  int<lower=1> nG; //number of regions
  int<lower=1> nS; //number of seasons
  int x[n]; //week of observation indicator
  int<lower = 1, upper = nG> group[n]; //indicator for each group
  int<lower = 0> y[n]; // observed value of flu
  int<lower = 1, upper = nS> seas[n]; // season indicator 
  int<lower=1> k; // number of mean parameters
  vector[k] mu0;
  cov_matrix[k] C0;
}

parameters {
  matrix[nS, k] ctheta [nG]; // Parameters for Each Year.  This is nG x nS x K array
}

transformed parameters {
  real lpsi[n]; // Logit of psi
  
  for(i in 1:n){
    if (x[i] < ctheta[group[i], seas[i], 4]){
      lpsi[i] = (ctheta[group[i], seas[i], 1] + (ctheta[group[i], seas[i], 3] - ctheta[group[i], seas[i], 1])*exp(-((x[i] - ctheta[group[i], seas[i], 4])^2)/(2*(exp(ctheta[group[i], seas[i], 5]))^2)));
    }
    else{
      lpsi[i] = (ctheta[group[i], seas[i], 2] + (ctheta[group[i], seas[i], 3] - ctheta[group[i], seas[i], 2])*exp(-((x[i] - ctheta[group[i], seas[i], 4])^2)/(2*(exp(ctheta[group[i], seas[i], 6]))^2)));
    }
    
  }
}

model {
  
  
  for(g in 1:nG){
    for(s in 1:nS){
      row(ctheta[g], s) ~ multi_normal(mu0, C0); // c(beta1[g], beta2[g], leta[g], mu[g], lsigma1[g], lsigma2[g])
    }
  }
  
  y ~ binomial_logit(num, lpsi);
}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = binomial_logit_lpmf(y[i] | num[i], lpsi[i]);
  }
}

