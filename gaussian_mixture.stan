data {
  int Nobs;
  vector[Nobs] xobs;

  real sigma_obs;

  real xth;

  int NNobs_max;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0> Lambda;
  real<lower=0> Lambda0;

  vector<lower=0>[Nobs] xobs_true;

  positive_ordered[NNobs_max] xnobs_true;
  vector<lower=0,upper=xth>[NNobs_max] xnobs;
}

model {
  /* Priors */
  mu ~ normal(0, 10);
  sigma ~ normal(0, 10);
  Lambda ~ normal(100.0, 100.0);
  Lambda0 ~ normal(0, NNobs_max);

  /* Observed likelihood */
  xobs ~ lognormal(log(xobs_true), sigma_obs);
  xobs_true ~ lognormal(mu, sigma);
  target += Nobs*log(Lambda);

  /* Non observed likelihood. */
  xnobs_true ~ lognormal(mu, sigma);
  for (i in 1:NNobs_max) {
    target += log_sum_exp(log(Lambda0) - log(xth),
			  log(Lambda) + lognormal_lpdf(xnobs[i] | log(xnobs_true[i]), sigma_obs));
  }

  target += -Lambda - Lambda0;
}
