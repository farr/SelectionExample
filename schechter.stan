functions {
  real schechter_lpdf(vector L, real alpha, real Lstar) {
    return sum(alpha*log(L) - L/Lstar - (alpha+1)*log(Lstar) - lgamma(alpha+1));
  }
}

data {
  int Nobs;
  vector[Nobs] Lobs;

  int NNobs_max;

  real Funcert; /* Fractional flux uncertainty */
  real zmax; /* Maximum redshift. */
  real Fth; /* Threshold flux */
}

parameters {
  /* Expected number of sources out to zmax. */
  real<lower=0> Lambda;

  /* Counter for the non-physical sources */
  real<lower=0> Lambda0;

  /* Power law at low luminosity. */
  real<lower=-1> alpha;

  /* Turnover to exponential decay. */
  real<lower=0> Lstar;

  /* True luminosity inferred for the observed systems. */
  vector<lower=0>[Nobs] Ltrue;

  /* Next parameters refer to the un-observed systems. */

  /* True luminosity of the (possibly) un-observed systems; making
     this vector `positive_ordered` makes the un-observed systems
     distinguishable. */
  positive_ordered[NNobs_max] Ltrue_nobs;
  /* True redshift of (possibly) unobserved systems. */
  vector<lower=0,upper=zmax>[NNobs_max] ztrue_nobs;

  /* To be non-observed, we must have a flux smaller than the flux
     limit. */
  vector<lower=0,upper=Fth>[NNobs_max] flux_nobs;
}

model {
  /* Priors. */
  Lambda ~ normal(0.0, 200.0);
  Lambda0 ~ normal(0.0, NNobs_max);
  alpha ~ normal(0.0, 1.0);
  Lstar ~ normal(0.0, 2.0);

  /* Observed systems (note that P_det == 1 for these systems, since
     their flux is above the limit). */
  Ltrue ~ schechter(alpha, Lstar);
  target += Nobs*log(Lambda);

  Lobs ~ lognormal(log(Ltrue), Funcert);

  /* Non-observed systems are a mix of *physical* systems (counted by
     Lambda), whose flux follows from their luminosity and
     *non-physical* systems (counted by Lambda0) whose flux is not
     tied to their luminosity at all, and instead distributed flat up
     to the flux threshold. */
  Ltrue_nobs ~ schechter(alpha, Lstar);
  /* Redshifts are flat */
  for (i in 1:NNobs_max) {
    real log_ex_flux = log(Ltrue_nobs[i]) - log(4.0*pi()) - 2.0*log(ztrue_nobs[i]);
    target += log_sum_exp(log(Lambda0) - log(Fth),
			  log(Lambda) + lognormal_lpdf(flux_nobs[i] | log_ex_flux, Funcert));
  }

  /* Poisson normalisation */
  target += -Lambda - Lambda0;
}
