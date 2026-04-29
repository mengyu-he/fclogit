// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ============================================================
// Utility functions
// ============================================================

// Inverse logit
inline double expit_scalar(const double eta) {
  if (eta >= 0.0) {
    const double z = std::exp(-eta);
    return 1.0 / (1.0 + z);
  } else {
    const double z = std::exp(eta);
    return z / (1.0 + z);
  }
}

vec expit_vec(const vec& eta) {
  vec out(eta.n_elem);
  for (int i = 0; i < eta.n_elem; ++i) out(i) = expit_scalar(eta(i));
  return out;
}

mat safe_inverse(const mat& A) {
  mat A_inv;
  bool ok = inv_sympd(A_inv, symmatu(A));
  if (!ok) A_inv = pinv(A);
  return A_inv;
}

// Convert a K*K*K cube to a K^2 * K matrix
// This makes the Firth term easy to evaluate for all k simultaneously
// [[Rcpp::export]]
mat append_slices(cube C) {
  mat res(C.n_rows * C.n_cols, C.n_slices, fill::zeros);
  for (int k = 0; k < C.n_slices; ++k) {
    res.col(k) = vectorise(C.slice(k));
  }
  return res;
}

// Outer product a %o% b %o% c stored as a cube with slice index from c
// [[Rcpp::export]]
cube outer_cube(vec a, vec b, vec c) {
  int n = a.n_elem;
  cube res(n, n, n, fill::zeros);
  for (int k = 0; k < n; ++k) {
    res.slice(k) = (a * b.t()) * c(k);
  }
  return res;
}

// Cube whose k-th slice is a * B[, k]^T
// [[Rcpp::export]]
cube outer_vec_mat(vec a, mat B) {
  int n = a.n_elem;
  cube res(n, n, n, fill::zeros);
  for (int k = 0; k < n; ++k) {
    res.slice(k) = a * B.col(k).t();
  }
  return res;
}

// Transpose each slice of a cube
// [[Rcpp::export]]
cube Apply_convert(cube C) {
  cube res(C.n_rows, C.n_cols, C.n_slices, fill::zeros);
  for (int k = 0; k < C.n_slices; ++k) {
    res.slice(k) = C.slice(k).t();
  }
  return res;
}

// Cube whose k-th slice is B * a[k]
// [[Rcpp::export]]
cube outer_mat_vec(mat B, vec a) {
  int n = a.n_elem;
  cube res(n, n, n, fill::zeros);
  for (int k = 0; k < n; ++k) {
    res.slice(k) = B * a(k);
  }
  return res;
}

// For an N*K matrix A, return the K*K*N cube whose i-th slice is A[i,]^T A[i,]
// [[Rcpp::export]]
cube outer_mat_row(mat A) {
  cube res(A.n_cols, A.n_cols, A.n_rows, fill::zeros);
  for (int i = 0; i < A.n_rows; ++i) {
    res.slice(i) = A.row(i).t() * A.row(i);
  }
  return res;
}

// ============================================================
// Howard-Gail recursions
// ============================================================

// (hg) 1st derivative of the HG normalizing constant
//
// Input:
//   beta: K-vector of regression coefficients
//   x: n x K covariate matrix for one stratum
//   m1: number of cases in this stratum
//
// Output:
//   mu = B_p(m1, n) / B(m1, n), the conditional expectation of 
//            sum_{selected i} z_i under the conditional logistic model
// [[Rcpp::export]]
List howard_gail_1(const vec beta, const mat x, double m1) {
  const int n_obs = x.n_rows;
  const int n_var = x.n_cols;
  const int M = static_cast<int>(m1);
  
  if (M <= 0 || M >= n_obs) stop("m1 must satisfy 0 < m1 < n in howard_gail_1().");
  
  vec psi = exp(x * beta);
  
  mat B(M, 2, fill::zeros);
  field<vec> B_r(M, 2);
  for (int m = 0; m < M; ++m) {
    for (int s = 0; s < 2; ++s) B_r(m, s) = zeros<vec>(n_var);
  }
  
  int this_val = 1;
  int last_val = 0;
  
  // Diagonal recursion: n = n1 + m. Only the current and previous
  // diagonals are stored, which is sufficient for the HG recurrence.
  for (int n1 = 0; n1 <= n_obs - M; ++n1) {
    std::swap(this_val, last_val);
    
    int m = 1;
    int n = n1 + m;
    rowvec z_i = x.row(n - 1);
    
    if (n == 1) {
      B(m - 1, this_val) = psi(n - 1);
      B_r(m - 1, this_val) = psi(n - 1) * z_i.t();
    } else {
      B(m - 1, this_val) = B(m - 1, last_val) + psi(n - 1);
      B_r(m - 1, this_val) = B_r(m - 1, last_val) + psi(n - 1) * z_i.t();
    }
    
    for (m = 2; m <= M; ++m) {
      n = n1 + m;
      z_i = x.row(n - 1);
      
      B(m - 1, this_val) = B(m - 1, last_val) + psi(n - 1) * B(m - 2, this_val);
      B_r(m - 1, this_val) = B_r(m - 1, last_val) +
        psi(n - 1) * B_r(m - 2, this_val) +
        psi(n - 1) * B(m - 2, this_val) * z_i.t();
    }
  }
  
  vec mu = B_r(M - 1, this_val) / B(M - 1, this_val);
  return List::create(Named("mu") = mu);
}

// (hg) recursions up to the 3rd derivative of the HG normalizing constant.
//
// Output:
//   mu: B_p / B.
//   sigma: B_pq / B - mu_p mu_q, i.e. the exact information contribution.
//   sigma_prime: derivative of sigma with respect to beta_r, stored slice-wise.
//                Slice r equals d sigma / d beta_r.
// [[Rcpp::export]]
List howard_gail_3(const vec beta, const mat x, double m1) {
  const int n_obs = x.n_rows;
  const int n_var = x.n_cols;
  const int M = static_cast<int>(m1);
  
  if (M <= 0 || M >= n_obs) stop("m1 must satisfy 0 < m1 < n in howard_gail_3().");
  
  vec psi = exp(x * beta);
  
  mat B(M, 2, fill::zeros);
  field<vec> B_r(M, 2);
  field<mat> B_rs(M, 2);
  field<cube> B_rst(M, 2);
  
  for (int m = 0; m < M; ++m) {
    for (int s = 0; s < 2; ++s) {
      B_r(m, s) = zeros<vec>(n_var);
      B_rs(m, s) = zeros<mat>(n_var, n_var);
      B_rst(m, s) = zeros<cube>(n_var, n_var, n_var);
    }
  }
  
  int this_val = 1;
  int last_val = 0;
  
  for (int n1 = 0; n1 <= n_obs - M; ++n1) {
    std::swap(this_val, last_val);
    
    int m = 1;
    int n = n1 + m;
    rowvec z_i = x.row(n - 1);
    vec z = z_i.t();
    cube z3 = outer_cube(z, z, z);
    
    if (n == 1) {
      B(m - 1, this_val) = psi(n - 1);
      B_r(m - 1, this_val) = psi(n - 1) * z;
      B_rs(m - 1, this_val) = psi(n - 1) * z * z.t();
      B_rst(m - 1, this_val) = psi(n - 1) * z3;
    } else {
      B(m - 1, this_val) = B(m - 1, last_val) + psi(n - 1);
      B_r(m - 1, this_val) = B_r(m - 1, last_val) + psi(n - 1) * z;
      B_rs(m - 1, this_val) = B_rs(m - 1, last_val) + psi(n - 1) * z * z.t();
      B_rst(m - 1, this_val) = B_rst(m - 1, last_val) + psi(n - 1) * z3;
    }
    
    for (m = 2; m <= M; ++m) {
      n = n1 + m;
      z_i = x.row(n - 1);
      z = z_i.t();
      z3 = outer_cube(z, z, z);
      
      B(m - 1, this_val) = B(m - 1, last_val) + psi(n - 1) * B(m - 2, this_val);
      
      B_r(m - 1, this_val) = B_r(m - 1, last_val) +
        psi(n - 1) * B_r(m - 2, this_val) +
        psi(n - 1) * B(m - 2, this_val) * z;
      
      B_rs(m - 1, this_val) = B_rs(m - 1, last_val) +
        psi(n - 1) * B_rs(m - 2, this_val) +
        psi(n - 1) * B(m - 2, this_val) * (z * z.t()) +
        psi(n - 1) * (z * B_r(m - 2, this_val).t() + B_r(m - 2, this_val) * z.t());
      
      cube tmp = outer_vec_mat(z, B_rs(m - 2, this_val));
      B_rst(m - 1, this_val) = B_rst(m - 1, last_val) +
        psi(n - 1) * B_rst(m - 2, this_val) +
        psi(n - 1) * B(m - 2, this_val) * z3 +
        psi(n - 1) * (outer_cube(z, z, B_r(m - 2, this_val)) +
        outer_cube(z, B_r(m - 2, this_val), z) +
        outer_cube(B_r(m - 2, this_val), z, z)) +
        psi(n - 1) * (tmp + Apply_convert(tmp) + outer_mat_vec(B_rs(m - 2, this_val), z));
    }
  }
  
  vec mu = B_r(M - 1, this_val) / B(M - 1, this_val);
  mat sigma = B_rs(M - 1, this_val) / B(M - 1, this_val) - mu * mu.t();
  cube sigma_prime = B_rst(M - 1, this_val) / B(M - 1, this_val) -
    outer_mat_vec(sigma, mu) -
    outer_vec_mat(mu, sigma) -
    Apply_convert(outer_vec_mat(mu, sigma)) -
    outer_cube(mu, mu, mu);
  
  return List::create(Named("mu") = mu,
                      Named("sigma") = sigma,
                      Named("sigma_prime") = sigma_prime);
}

// (hg) score over all strata.
//
// Note:
//   x_all is ordered by stratum, and within each stratum the cases appear
//   before controls. This is the ordering created by the R function.
// [[Rcpp::export]]
vec score_hg(vec beta, mat x_all, vec ns, vec m1s, bool firth) {
  const int n_strata = ns.n_elem;
  const int n_var = x_all.n_cols;
  vec score(n_var, fill::zeros);
  
  mat I(n_var, n_var, fill::zeros);
  cube I_prime(n_var, n_var, n_var, fill::zeros);
  
  int low = 0;
  for (int s = 0; s < n_strata; ++s) {
    const int n_s = static_cast<int>(ns(s));
    const int m1_s = static_cast<int>(m1s(s));
    const int up = low + n_s - 1;
    
    mat x_s = x_all.rows(low, up);
    vec y_s = join_cols(ones<vec>(m1_s), zeros<vec>(n_s - m1_s));
    
    if (!firth) {
      List res = howard_gail_1(beta, x_s, m1_s);
      vec mu = res["mu"];
      score += x_s.t() * y_s - mu;
    } else {
      List res = howard_gail_3(beta, x_s, m1_s);
      vec mu = res["mu"];
      mat I_s = res["sigma"];
      cube I_prime_s = res["sigma_prime"];
      
      score += x_s.t() * y_s - mu;
      I += I_s;
      I_prime += I_prime_s;
    }
    
    low = up + 1;
  }
  
  if (firth) {
    rowvec I_inv_vec = vectorise(safe_inverse(I), 1);
    mat I_prime_mat = append_slices(I_prime);
    score += trans(0.5 * I_inv_vec * I_prime_mat);
  }
  
  return score;
}

// ============================================================
// Saddlepoint approximation
// ============================================================

// (sp) Alpha equation for a grouped stratum: m1 - sum_j c_j p_j = 0.
// [[Rcpp::export]]
double eqn_alpha_sp(vec p, vec c, double m1) {
  return m1 - dot(c, p);
}

// (sp) SP beta-score contribution for one grouped stratum.
vec sp_score_beta_one_stratum(const vec& beta,
                              const double alpha,
                              const mat& x_unique_j,
                              const vec& c,
                              const vec& t1,
                              const double n) {
  vec eta = alpha + x_unique_j * beta;
  vec p = expit_vec(eta);
  vec w = c % p % (1.0 - p);
  double sum_w = sum(w);
  if (sum_w <= 0.0 || !std::isfinite(sum_w)) stop("SP failed: non-positive or non-finite sum of weights.");
  
  rowvec x_bar = (w.t() * x_unique_j) / sum_w;
  mat z_minus_z_bar = x_unique_j.each_row() - x_bar;
  
  // Equivalent to - n/(n-1) * sum_j (w_j/sum_w) * (p_j - p_bar) * (z_j-z_bar),
  // because sum_j w_j (z_j-z_bar) = 0.
  vec correction = (n / (n - 1.0) / sum_w) * (z_minus_z_bar.t() * (w % p));
  return x_unique_j.t() * (t1 - c % p) - correction;
}

// (sp) current analytic SP information contribution
mat sp_info_one_stratum_current(const vec& beta,
                                const double alpha,
                                const mat& x_unique_j,
                                const vec& c,
                                const double n) {
  const int n_var = x_unique_j.n_cols;
  
  vec eta = alpha + x_unique_j * beta;
  vec p = expit_vec(eta);
  vec w = c % p % (1.0 - p);
  double sum_w = sum(w);
  if (sum_w <= 0.0 || !std::isfinite(sum_w)) stop("SP failed: non-positive or non-finite sum of weights.");
  
  rowvec x_bar = (w.t() * x_unique_j) / sum_w;
  mat z_minus_z_bar = x_unique_j.each_row() - x_bar;
  double p_bar = dot(w, p) / sum_w;
  
  mat s1(n_var, n_var, fill::zeros);
  mat s2(n_var, n_var, fill::zeros);
  vec s3(n_var, fill::zeros);
  cube z_outer = outer_mat_row(z_minus_z_bar);
  
  vec coef1 = w % ((1.0 - 2.0 * p) % (p - p_bar) + w / c);
  vec coef2 = w % (1.0 - 2.0 * p);
  
  for (int k = 0; k < x_unique_j.n_rows; ++k) {
    s1 += w(k) * z_outer.slice(k);
    s2 += coef1(k) * z_outer.slice(k);
    s3 += coef2(k) * z_minus_z_bar.row(k).t();
  }
  
  return s1 + (n / (n - 1.0) / sum_w) * s2 +
    (n / (n - 1.0) / (2.0 * std::pow(sum_w, 2.0))) * s3 * s3.t();
}

// (sp) score function. 
// Note:
//   The unknown vector betas contains all stratum intercepts followed by the common beta vector:
//   x_unique, cs, and t1s are grouped by stratum. 
// [[Rcpp::export]]
vec score_sp(vec betas, mat x_unique, vec ns, vec m1s, vec t1s, vec cs,
             vec first_unique, vec last_unique, bool firth) {
  const int n_var = x_unique.n_cols;
  const int n_strata = ns.n_elem;
  
  vec alphas = betas.subvec(0, n_strata - 1);
  vec beta = betas.subvec(n_strata, n_strata + n_var - 1);
  vec score(n_var + n_strata, fill::zeros);
  
  mat I(n_var, n_var, fill::zeros);
  cube I_prime(n_var, n_var, n_var, fill::zeros);
  
  for (int s = 0; s < n_strata; ++s) {
    const int ind_first = static_cast<int>(first_unique(s)) - 1;
    const int ind_last  = static_cast<int>(last_unique(s)) - 1;
    const double n_s = ns(s);
    const double m1_s = m1s(s);
    const double alpha_s = alphas(s);
    
    mat x_s = x_unique.rows(ind_first, ind_last);
    vec c_s = cs.subvec(ind_first, ind_last);
    vec t1_s = t1s.subvec(ind_first, ind_last);
    
    vec eta = alpha_s + x_s * beta;
    vec p = expit_vec(eta);
    
    score(s) += eqn_alpha_sp(p, c_s, m1_s);
    score.subvec(n_strata, n_strata + n_var - 1) +=
      sp_score_beta_one_stratum(beta, alpha_s, x_s, c_s, t1_s, n_s);
    
    if (firth) {
      I += sp_info_one_stratum_current(beta, alpha_s, x_s, c_s, n_s);
      
      // Analytic derivative of the current SP information expression
      const int J_s = x_s.n_rows;
      vec w = c_s % p % (1.0 - p);
      double sum_w = sum(w);
      rowvec x_bar = (w.t() * x_s) / sum_w;
      mat z_minus_z_bar = x_s.each_row() - x_bar;
      double p_bar = dot(w, p) / sum_w;
      cube z_outer = outer_mat_row(z_minus_z_bar);
      
      vec coef1 = w % ((1.0 - 2.0 * p) % (p - p_bar) + w / c_s);
      vec coef2 = w % (1.0 - 2.0 * p);
      
      mat s2(n_var, n_var, fill::zeros);
      vec s3(n_var, fill::zeros);
      for (int j = 0; j < J_s; ++j) {
        s2 += coef1(j) * z_outer.slice(j);
        s3 += coef2(j) * z_minus_z_bar.row(j).t();
      }
      
      mat zbar_prime = z_minus_z_bar.t() * (z_minus_z_bar.each_col() % coef2) / sum_w;
      vec pbar_prime = z_minus_z_bar.t() * coef1 / sum_w;
      
      for (int kpp = 0; kpp < n_var; ++kpp) {
        mat tt1(n_var, n_var, fill::zeros);
        mat t2(n_var, n_var, fill::zeros);
        mat t3(n_var, n_var, fill::zeros);
        
        for (int j = 0; j < J_s; ++j) {
          mat z_outer_j = z_outer.slice(j);
          mat weighted_z_outer = coef2(j) * z_minus_z_bar(j, kpp) * z_outer_j;
          
          tt1 += weighted_z_outer;
          t2 += ((1.0 - 2.0 * p(j)) * (p(j) - p_bar) + 2.0 * w(j) / c_s(j)) * weighted_z_outer +
            w(j) * (p(j) * (1.0 - p(j)) * z_minus_z_bar(j, kpp) *
            (1.0 - 4.0 * p(j) + 2.0 * p_bar) -
            (1.0 - 2.0 * p(j)) * pbar_prime(kpp)) * z_outer_j +
            w(j) * ((1.0 - 2.0 * p(j)) * (p(j) - p_bar) + w(j) / c_s(j)) *
            (-z_minus_z_bar.row(j).t() * zbar_prime.col(kpp).t() -
            zbar_prime.col(kpp) * z_minus_z_bar.row(j));
        }
        
        t2 = (n_s / (n_s - 1.0)) * (t2 / sum_w - s2 * s3(kpp) / std::pow(sum_w, 2.0));
        t3 = (n_s / (n_s - 1.0) / 2.0) *
          (-2.0 * s2.row(kpp).t() / sum_w - s3(kpp) / std::pow(sum_w, 2.0)) *
          s3.t() / sum_w;
        t3 = t3 + t3.t();
        
        I_prime.slice(kpp) += tt1 + t2 + t3;
      }
    }
  }
  
  if (firth) {
    rowvec I_inv_vec = vectorise(safe_inverse(I), 1);
    mat I_prime_mat = append_slices(I_prime);
    score.subvec(n_strata, n_strata + n_var - 1) +=
      trans(0.5 * I_inv_vec * I_prime_mat);
  }
  
  return score;
}

// (sp) inverse Fisher information based on the current analytic SP information.
// [[Rcpp::export]]
mat get_fisher_inv(vec betas, mat x_unique, vec ns, vec m1s, vec t1s, vec cs,
                   vec first_unique, vec last_unique) {
  const int n_var = x_unique.n_cols;
  const int n_strata = ns.n_elem;
  vec alphas = betas.subvec(0, n_strata - 1);
  vec beta = betas.subvec(n_strata, n_strata + n_var - 1);
  
  mat I(n_var, n_var, fill::zeros);
  for (int s = 0; s < n_strata; ++s) {
    const int ind_first = static_cast<int>(first_unique(s)) - 1;
    const int ind_last  = static_cast<int>(last_unique(s)) - 1;
    mat x_s = x_unique.rows(ind_first, ind_last);
    vec c_s = cs.subvec(ind_first, ind_last);
    I += sp_info_one_stratum_current(beta, alphas(s), x_s, c_s, ns(s));
  }
  
  return safe_inverse(I);
}

// ============================================================
// Hybrid HG/SP score used by the R wrapper
// ============================================================

// stratum_method: 1 = HG, 2 = SP, one entry per stratum.
// hg_flip: 1 if the HG contribution is computed after switching cases and 
// controls within that stratum.  In that case the HG recursion is evaluated 
// at -beta with m0 cases, and the beta-score and dI/dbeta are multiplied by -1
// when converted back to the original parameterization.
//
// [[Rcpp::export]]
vec score_hybrid(vec par, mat x_all, mat x_unique, vec ns, vec m1s,
                 vec t1s, vec cs, vec first_unique, vec last_unique,
                 vec stratum_method, vec hg_flip, bool firth) {
  const int n_strata = ns.n_elem;
  const int n_var = x_all.n_cols;
  int n_sp = 0;
  for (int ss = 0; ss < n_strata; ++ss) {
    if (static_cast<int>(stratum_method(ss)) == 2) n_sp++;
  }
  
  vec alphas;
  if (n_sp > 0) alphas = par.subvec(0, n_sp - 1);
  vec beta = par.subvec(n_sp, n_sp + n_var - 1);
  
  vec score(n_sp + n_var, fill::zeros);
  mat I(n_var, n_var, fill::zeros);
  cube I_prime(n_var, n_var, n_var, fill::zeros);
  
  int low = 0;
  int sp_ind = 0;
  
  for (int s = 0; s < n_strata; ++s) {
    const int n_s = static_cast<int>(ns(s));
    const int m1_s = static_cast<int>(m1s(s));
    const int up = low + n_s - 1;
    
    if (static_cast<int>(stratum_method(s)) == 1) {
      mat x_s = x_all.rows(low, up);
      bool flip_s = (static_cast<int>(hg_flip(s)) == 1);
      int M = flip_s ? (n_s - m1_s) : m1_s;
      double sign_s = flip_s ? -1.0 : 1.0;
      vec beta_hg = flip_s ? -beta : beta;
      
      vec x_y(n_var, fill::zeros);
      if (!flip_s) {
        if (m1_s > 0) x_y = trans(x_s.rows(0, m1_s - 1)) * ones<vec>(m1_s);
      } else {
        if (n_s - m1_s > 0) x_y = trans(x_s.rows(m1_s, n_s - 1)) * ones<vec>(n_s - m1_s);
      }
      
      if (!firth) {
        List res = howard_gail_1(beta_hg, x_s, M);
        vec mu = res["mu"];
        score.subvec(n_sp, n_sp + n_var - 1) += sign_s * (x_y - mu);
      } else {
        List res = howard_gail_3(beta_hg, x_s, M);
        vec mu = res["mu"];
        mat I_s = res["sigma"];
        cube I_prime_s = res["sigma_prime"];
        
        score.subvec(n_sp, n_sp + n_var - 1) += sign_s * (x_y - mu);
        I += I_s;
        I_prime += sign_s * I_prime_s;
      }
      
    } else {
      const int ind_first = static_cast<int>(first_unique(s)) - 1;
      const int ind_last  = static_cast<int>(last_unique(s)) - 1;
      const double n_sd = ns(s);
      const double m1_sd = m1s(s);
      const double alpha_s = alphas(sp_ind);
      
      mat x_s = x_unique.rows(ind_first, ind_last);
      vec c_s = cs.subvec(ind_first, ind_last);
      vec t1_s = t1s.subvec(ind_first, ind_last);
      
      vec eta = alpha_s + x_s * beta;
      vec p = expit_vec(eta);
      
      score(sp_ind) = eqn_alpha_sp(p, c_s, m1_sd);
      score.subvec(n_sp, n_sp + n_var - 1) +=
        sp_score_beta_one_stratum(beta, alpha_s, x_s, c_s, t1_s, n_sd);
      
      if (firth) {
        I += sp_info_one_stratum_current(beta, alpha_s, x_s, c_s, n_sd);
        
        const int J_s = x_s.n_rows;
        vec w = c_s % p % (1.0 - p);
        double sum_w = sum(w);
        rowvec x_bar = (w.t() * x_s) / sum_w;
        mat z_minus_z_bar = x_s.each_row() - x_bar;
        double p_bar = dot(w, p) / sum_w;
        cube z_outer = outer_mat_row(z_minus_z_bar);
        
        vec coef1 = w % ((1.0 - 2.0 * p) % (p - p_bar) + w / c_s);
        vec coef2 = w % (1.0 - 2.0 * p);
        
        mat s2(n_var, n_var, fill::zeros);
        vec s3(n_var, fill::zeros);
        for (int j = 0; j < J_s; ++j) {
          s2 += coef1(j) * z_outer.slice(j);
          s3 += coef2(j) * z_minus_z_bar.row(j).t();
        }
        
        mat zbar_prime = z_minus_z_bar.t() * (z_minus_z_bar.each_col() % coef2) / sum_w;
        vec pbar_prime = z_minus_z_bar.t() * coef1 / sum_w;
        
        for (int kpp = 0; kpp < n_var; ++kpp) {
          mat tt1(n_var, n_var, fill::zeros);
          mat t2(n_var, n_var, fill::zeros);
          mat t3(n_var, n_var, fill::zeros);
          
          for (int j = 0; j < J_s; ++j) {
            mat z_outer_j = z_outer.slice(j);
            mat weighted_z_outer = coef2(j) * z_minus_z_bar(j, kpp) * z_outer_j;
            
            tt1 += weighted_z_outer;
            t2 += ((1.0 - 2.0 * p(j)) * (p(j) - p_bar) + 2.0 * w(j) / c_s(j)) * weighted_z_outer +
              w(j) * (p(j) * (1.0 - p(j)) * z_minus_z_bar(j, kpp) *
              (1.0 - 4.0 * p(j) + 2.0 * p_bar) -
              (1.0 - 2.0 * p(j)) * pbar_prime(kpp)) * z_outer_j +
              w(j) * ((1.0 - 2.0 * p(j)) * (p(j) - p_bar) + w(j) / c_s(j)) *
              (-z_minus_z_bar.row(j).t() * zbar_prime.col(kpp).t() -
              zbar_prime.col(kpp) * z_minus_z_bar.row(j));
          }
          
          t2 = (n_sd / (n_sd - 1.0)) * (t2 / sum_w - s2 * s3(kpp) / std::pow(sum_w, 2.0));
          t3 = (n_sd / (n_sd - 1.0) / 2.0) *
            (-2.0 * s2.row(kpp).t() / sum_w - s3(kpp) / std::pow(sum_w, 2.0)) *
            s3.t() / sum_w;
          t3 = t3 + t3.t();
          
          I_prime.slice(kpp) += tt1 + t2 + t3;
        }
      }
      
      sp_ind++;
    }
    
    low = up + 1;
  }
  
  if (firth) {
    rowvec I_inv_vec = vectorise(safe_inverse(I), 1);
    mat I_prime_mat = append_slices(I_prime);
    score.subvec(n_sp, n_sp + n_var - 1) +=
      trans(0.5 * I_inv_vec * I_prime_mat);
  }
  
  return score;
}
