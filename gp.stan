#include covs.stan // definitions for the covariance functions
data {
  int<lower=1> N1;
  int<lower=1> N2;
  int<lower=1> D_w;
  vector[N1] x1_d;          // days (training set)
  vector[D_w] x1_w[N1];     // weather-data (training set)
  vector[N1] y1;            // cycling count (training set)
  vector[N2] x2_d;          // days (all)
  vector[D_w] x2_w[N2];     // weather-data (all)
  real<lower=0> delta;      // regularization
  real<lower=0> sd_sig;     // sd for sigma
  real<lower=0> sd_a;       // sd for marginal standard deviations
  real<lower=0> gamma;      // period of the periodic CF
}
parameters {
  real<lower=0> alpha_pe;       // marginal sd for the periodic component
  real<lower=0> alpha_ard;      // marginal sd for the weather-data component
  real<lower=0> rho_pe;         // length-scale for the periodic component
  vector<lower=0>[D_w] rho_ard; // length-scales for the weather-data component
  real<lower=0> sigma;          // noise in the observations
  vector[N1] eta;               // for the cholesky transformed implementation
}
model {
  vector[N1] f1;
  {
    matrix[N1, N1] K_pe = cov_periodic_s(x1_d, alpha_pe, rho_pe, gamma);
    matrix[N1, N1] K_ard = cov_ard_s(x1_w, alpha_ard, rho_ard);
    matrix[N1, N1] K = K_pe + K_ard + diag_matrix(rep_vector(delta, N1));
    matrix[N1, N1] L_K = cholesky_decompose(K);

    f1 = L_K * eta;
  }

  alpha_ard ~ normal(0, sd_a);
  alpha_pe ~ normal(0, sd_a);
  rho_ard ~ inv_gamma(5, 5);
  rho_pe ~ inv_gamma(5, 5);
  sigma ~ normal(0, sd_sig);
  eta ~ normal(0, 1);
  y1 ~ normal(f1, sigma);
}
generated quantities {
  vector[N2] f2;        // latent function at all points
  vector[N2] y2;        // predictions of y at all points
  {
    matrix[N1, N1] K11 = cov_periodic_s(x1_d, alpha_pe, rho_pe, gamma) +
                          cov_ard_s(x1_w, alpha_ard, rho_ard);
    matrix[N1, N2] K12 = cov_periodic(x1_d, x2_d, alpha_pe, rho_pe, gamma) +
                          cov_ard(x1_w, x2_w, alpha_ard, rho_ard);
    matrix[N2, N2] K22 = cov_periodic_s(x2_d, alpha_pe, rho_pe, gamma) +
                          cov_ard_s(x2_w, alpha_ard, rho_ard);

    f2 = gp_pred_rng(y1, K11, K12, K22, sigma, delta);
  }

  for (i in 1:N2) y2[i] = normal_rng(f2[i], sigma);

}
