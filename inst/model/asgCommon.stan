data {
  int<lower=1> n; //total number of obs over all locations and seasons
  int<lower=1> num[n]; //count of total patients
  int x[n]; //week of observation indicator
  int<lower = 0> y[n]; // observed value of flu
}

parameters {
  real beta1;
  real beta2;
  real nu;
  real mu;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
}

transformed parameters {
  real lpsi[n]; // Logit of psi
  
  for(i in 1:n){
    if (x[i] < mu){
      lpsi[i] = (beta1 + (nu - beta1)*exp(-((x[i] - mu)^2)/(2*(sigma1)^2)));
    }
    else{
      lpsi[i] = (beta2 + (nu - beta2)*exp(-((x[i] - mu)^2)/(2*(sigma2)^2)));
    }
    
  }
}

model {
  beta1 ~ normal(0,10);
  beta2 ~ normal(0,10);
  nu ~ normal(0, 5);
  mu ~ normal(15, 4);
  sigma1 ~ student_t(4, 0, 1);
  sigma2 ~ student_t(4, 0, 1);
  
  y ~ binomial_logit(num, lpsi);
}

generated quantities {
  vector[n] log_lik;
  for(i in 1:n){
    log_lik[i] = binomial_logit_lpmf(y[i] | num[i], lpsi[i]);
  }
}

