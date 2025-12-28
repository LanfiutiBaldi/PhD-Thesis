// generated with brms 2.19.0
functions {
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> K_sigma;  // number of population-level effects
  matrix[N, K_sigma] X_sigma;  // population-level design matrix
  int<lower=1> K_alpha;  // number of population-level effects
  matrix[N, K_alpha] X_alpha;  // population-level design matrix
  int prior_only;  // should the likelihood be ignored?
  int<lower=1> nage; // number of ages
  int<lower=1> nperiod; // number of periods 
  int<lower=1> ncohort; // number of cohorts 
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b_raw;  // population-level effects
  real Intercept;  // temporary intercept for centered predictors
  vector[K_sigma] b_sigma;  // population-level effects
  vector[K_alpha] b_alpha;  // population-level effects
}
transformed parameters {
  vector[nage-1] b_1 = b_raw[1:(nage-1)];  
  vector[nperiod-1] b_2 = b_raw[nage:(nage+nperiod-2)];  
  vector[ncohort-1] b_3 = b_raw[(nage+nperiod-1):(nage+nperiod+ncohort-3)];
  vector[Kc] b;
  b = append_row(b_1, append_row(b_2, b_3));
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0.6, 2.5);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] sigma = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] alpha = rep_vector(0.0, N);
    // parameters used to transform the skew-normal distribution
    vector[N] delta;  // transformed alpha parameter
    vector[N] omega;  // scale parameter
    mu += Intercept + Xc * b;
    sigma += X_sigma * b_sigma;
    alpha += X_alpha * b_alpha;
    sigma = exp(sigma);
    // use efficient skew-normal parameterization
    for (n in 1:N) {
      delta[n] = alpha[n] / sqrt(1 + alpha[n]^2);
      omega[n] = sigma[n] / sqrt(1 - sqrt(2 / pi())^2 * delta[n]^2);
      mu[n] = mu[n] - omega[n] * delta[n] * sqrt(2 / pi());
    }
    target += skew_normal_lpdf(Y | mu, omega, alpha);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
}
