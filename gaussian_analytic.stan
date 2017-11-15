data {
  int Nobs;
  vector[Nobs] xobs;

  real sigma_obs;

  real xth;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> Lambda;

  vector<lower=0>[Nobs] xobs_true;
}

model {
  /* Priors */
  mu ~ normal(0, 10);
  sigma ~ normal(0, 10);
  Lambda ~ normal(100.0, 100.0);

  /* Observed likelihood */
  xobs ~ lognormal(log(xobs_true), sigma_obs);
  xobs_true ~ lognormal(mu, sigma);
  target += Nobs*log(Lambda);

  target += -Lambda*exp(lognormal_lccdf(xth | mu, sqrt(sigma_obs*sigma_obs + sigma*sigma)));
}
