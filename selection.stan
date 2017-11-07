data {
  int Nobs;
  vector[Nobs] Lobs;

  int Nmax; /* Maximum number of unobserved sources. */

  real zmax; /* Maximum redshift. */
  real Fth; /* Threshold flux. */
  real Funcert; /* Fractional flux uncertainty */
}

parameters {
  real<lower=0> N;
  real<lower=-1> alpha;
  real<lower=0> Lstar;
  vector<lower=0>[Nobs] Ltrue;
  vector<lower=0>[Nmax] Ltrue_unobs;
  vector<lower=0,upper=zmax>[Nmax] ztrue_unobs; 
}

model {
  vector[Nmax] ex_logflux;
  vector[Nmax] log_Pndet;
  real log_p_unobs_acc_term;
  real log_p_unobs_acc;

  /* Priors on population parameters (I made these up to be reasonable). */
  alpha ~ normal(0, 1);
  Lstar ~ normal(0, 1);
  N ~ normal(0, (Nobs+Nmax)/3.0); /* Nmax is a 3-sigma upward fluctuation by prior. */
  
  /* Schechter distribution for the (unobserved) true Luminosities. */
  target += sum(log(N) + alpha*log(Ltrue) - Ltrue/Lstar - (alpha+1)*log(Lstar) - lgamma(1+alpha));

  /* Likelihood for the observed lumonisity. */
  Lobs ~ lognormal(log(Ltrue), Funcert);

  /* Now for the non-detections. */

  /* Redshifts are distributed flat, so no term needed for
     z-distribution. */
    
  /* The probability of non-detection is just the fraction of the
     log-normal that lies below the flux threshold. */
  ex_logflux = log(Ltrue_unobs) - log(4.0*pi()) - 2.0*log(ztrue_unobs);
  for (i in 1:Nmax) {
    log_Pndet[i] = normal_lcdf(log(Fth) | ex_logflux[i], Funcert);
  }

  /* Schechter term for each unobserved galaxy. */
  target += sum(alpha*log(Ltrue_unobs) - Ltrue_unobs/Lstar - (alpha+1)*log(Lstar) - lgamma(1+alpha));

  /* Now we sum of the (unknown) number of non-detections: 

     log_p = log(1.0 + <one non-detection> + <two non-detections> + ...)
  */
  log_p_unobs_acc = 0.0;
  log_p_unobs_acc_term = 0.0;
  for (i in 1:Nmax) {
    /* <N+1 non-detections> = <N non-detections>*<N+1 term> */
    log_p_unobs_acc_term = log_p_unobs_acc_term + log_Pndet[i] + log(N) - log(i);

    log_p_unobs_acc = log_sum_exp(log_p_unobs_acc, log_p_unobs_acc_term);
  }
  target += log_p_unobs_acc;

  target += -N; /* exp(-N) normalisation for Poisson. */   
}
