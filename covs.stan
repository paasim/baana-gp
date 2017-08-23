functions {

  real k_sexp(real x1, real x2, real rho) {
    return exp(-0.5*square((x1 - x2) / rho));
  }
  matrix cov_sexp_s(vector x, real alpha, real rho) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_sexp(x[i], x[j], rho);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  matrix cov_sexp(vector x1, vector x2, real alpha, real rho) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_sexp(x1[i], x2[j], rho);

    return K;
  }

  real k_ard(vector x1, vector x2, vector rho) {
    return exp(-0.5*dot_self((x1 - x2) ./ rho));
  }
  matrix cov_ard_s(vector[] x, real alpha, vector rho) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_ard(x[i], x[j], rho);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  matrix cov_ard(vector[] x1, vector[] x2, real alpha, vector rho) {
    int N1 = size(x1);
    int N2 = size(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_ard(x1[i], x2[j], rho);

    return K;
  }

  real k_exp(real x1, real x2, real rho) {
    return exp(-fabs(x1 - x2) / rho);
  }
  matrix cov_exp_s(vector x, real alpha, real rho) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_exp(x[i], x[j], rho);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  matrix cov_exp(vector x1, vector x2, real alpha, real rho) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_exp(x1[i], x2[j], rho);

    return K;
  }

  real k_ma32(real x1, real x2, real rho) {
    real sq3r_rho = sqrt(3)*fabs(x1 - x2)/rho;
    return (1+sq3r_rho) * exp(-sq3r_rho);
  }
  matrix cov_ma32_s(vector x, real alpha, real rho) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_ma32(x[i], x[j], rho);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  matrix cov_ma32(vector x1, vector x2, real alpha, real rho) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_ma32(x1[i], x2[j], rho);

    return K;
  }

  real k_ma52(real x1, real x2, real rho) {
    real sq5r_rho = sqrt(5)*fabs(x1 - x2)/rho;
    return (1+sq5r_rho +square(sq5r_rho)/3) * exp(-sq5r_rho);
  }
  matrix cov_ma52_s(vector x, real alpha, real rho) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_ma52(x[i], x[j], rho);
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
  matrix cov_ma52(vector x1, vector x2, real alpha, real rho) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_ma52(x1[i], x2[j], rho);

    return K;
  }

  real k_periodic(real x1, real x2, real rho, real gamma) {
    return exp(-2*square(sin(pi()*(x1 - x2)/gamma) / rho));
  }
  matrix cov_periodic_s(vector x, real alpha, real rho, real gamma) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i,i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_periodic(x[i], x[j], rho, gamma);
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = sq_alpha;
    return K;
  }
  matrix cov_periodic(vector x1, vector x2, real alpha, real rho, real gamma) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_periodic(x1[i], x2[j], rho, gamma);

    return K;
  }

  real k_lperiodic(real x1, real x2, real rho_se, real rho_pe, real gamma) {
    return exp(-0.5*square((x1 - x2) / rho_se)
                -2*square(sin(pi()*(x1 - x2)/gamma) / rho_pe));
  }
  matrix cov_lperiodic_s(vector x, real alpha, real rho_se, real rho_pe, real gamma) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i,i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_lperiodic(x[i], x[j], rho_se, rho_pe, gamma);
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = sq_alpha;
    return K;
  }
  matrix cov_lperiodic(vector x1, vector x2, real alpha, real rho_se, real rho_pe, real gamma) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_lperiodic(x1[i], x2[j], rho_se, rho_pe, gamma);

    return K;
  }

  real k_pol(real x1, real x2, real b, int p) {
    return pow(x1*x2 + b^2, p);
  }
  matrix cov_pol_s(vector x, real alpha, real b, int p) {
    int N = num_elements(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    for (i in 1:(N-1)) {
      K[i,i] = sq_alpha*k_pol(x[i], x[i], b, p);
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha*k_pol(x[i], x[j], b, p);
        K[j, i] = K[i, j];
      }
    }
    K[N,N] = sq_alpha*k_pol(x[N], x[N], b, p);
    return K;
  }
  matrix cov_pol(vector x1, vector x2, real alpha, real b, int p) {
    int N1 = num_elements(x1);
    int N2 = num_elements(x2);
    matrix[N1, N2] K;
    real sq_alpha = square(alpha);
    for (i in 1:N1)
      for (j in 1:N2)
        K[i, j] = sq_alpha*k_pol(x1[i], x2[j], b, p);

    return K;
  }

  vector gp_pred_rng(vector y1, matrix K11, matrix K12, matrix K22,
                     real sigma, real delta) {
    int N1 = rows(K11);
    int N2 = rows(K22);
    vector[N2] f2;
    {
      matrix[N1, N1] sq_sig = diag_matrix(rep_vector(square(sigma), N1));
      matrix[N1, N1] L_K = cholesky_decompose(K11 + sq_sig);
      matrix[N1, N2] v_pred = mdivide_left_tri_low(L_K, K12);
      matrix[N2, N2] cov_f2 = K22 - v_pred' * v_pred;
      vector[N1] K_div_y1 = mdivide_left_tri_low(L_K, y1);
      vector[N2] f2_mu = (K12' * mdivide_right_tri_low(K_div_y1', L_K)');
      f2 = multi_normal_rng(f2_mu, cov_f2 + diag_matrix(rep_vector(delta, N2)));
    }
    return f2;
  }
}
