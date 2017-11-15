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

  vector<lower=0>[Nobs] xobs_true;

  vector<lower=0>[NNobs_max] xnobs_true;
  vector<lower=0,upper=xth>[NNobs_max] xnobs;
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

  /* Non observed likelihood. */
  xnobs_true ~ lognormal(mu, sigma);
  target += -NNobs_max*log(xth); /* Default flat prior for xnobs */

  {
    vector[NNobs_max+1] log_poisson_term;

    for (i in 1:NNobs_max) {
      log_poisson_term[i+1] = log(Lambda) + lognormal_lpdf(xnobs[i] | log(xnobs_true[i]), sigma_obs) - log(i) + log(xth);
    }
    log_poisson_term[1] = 0.0;

    log_poisson_term = cumulative_sum(log_poisson_term);

    target += log_sum_exp(log_poisson_term);
  }

  target += -Lambda;
}
