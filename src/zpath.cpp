#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

inline double digamma_cpp(double x) { return R::digamma(x); }
inline double trigamma_cpp(double x) { return R::trigamma(x); }
inline double lgamma_cpp(double x) { return R::lgammafn(x); }

// [[Rcpp::export]]
Rcpp::List zpath(arma::mat beta,
                 double phi,
                 const arma::mat &logy,
                 const arma::mat &y,
                 const arma::mat &x,
                 const Rcpp::List &txi,
                 const Rcpp::List &Pidx,
                 double lambda,
                 int n,
                 int D,
                 int p,
                 int K,
                 const arma::vec &pen,
                 int maxit,
                 double tol) {
  
  const int dim = K * p;
  
  // lambda_vec and lamI exactly as in R
  arma::vec lambda_vec(dim);
  for (int j = 0; j < p; ++j)
    for (int k = 0; k < K; ++k)
      lambda_vec[k + K * j] = lambda * pen[k];
  
  arma::mat lamI = arma::diagmat(lambda_vec);
  
  // match R: phi reset inside .zpath
  phi = 1.0;
  double loglik_old = -std::numeric_limits<double>::infinity();
  double loglik_pen = NA_REAL;
  
  for (int iter = 0; iter < maxit; ++iter) {
    
    // eta = x %*% beta
    arma::mat eta = x * beta;
    arma::mat exp_eta = arma::exp(eta);
    
    // mu = cbind(1, exp_eta) / (1 + rowsums(exp_eta))
    arma::vec sum_exp = arma::sum(exp_eta, 1);
    arma::mat mu(n, D);
    mu.col(0) = 1.0 / (1.0 + sum_exp);
    for (int j = 0; j < p; ++j)
      mu.col(j + 1) = exp_eta.col(j) / (1.0 + sum_exp);
    
    double psi1_phi = trigamma_cpp(phi);
    double loglik = 0.0;
    arma::vec S_vec(dim, fill::zeros);
    arma::mat I_mat(dim, dim, fill::zeros);
    double S_phi = 0.0;
    double H_phi = 0.0;
    
    for (int i = 0; i < n; ++i) {
      
      // Pi indices (1-based in R, convert to 0-based)
      IntegerVector Pi_R = Pidx[i];
      int mPi = Pi_R.size();
      arma::uvec Pi(mPi);
      for (int idx = 0; idx < mPi; ++idx)
        Pi[idx] = Pi_R[idx] - 1;
      
      arma::rowvec x_i = x.row(i);
      arma::rowvec mu_i = mu.row(i);
      
      // Mi, mu2, alpha2, y_log
      double Mi = 0.0;
      for (int idx = 0; idx < mPi; ++idx)
        Mi += mu_i[Pi[idx]];
      
      arma::vec mu2(mPi);
      arma::vec alpha2(mPi);
      arma::vec y_log(mPi);
      for (int idx = 0; idx < mPi; ++idx) {
        double val = mu_i[Pi[idx]] / Mi;
        mu2[idx]    = val;
        alpha2[idx] = phi * val;
        y_log[idx]  = logy(i, Pi[idx]);
      }
      
      // loglik, S_phi, H_phi (Dirichlet part)
      double sum_lg_alpha2   = 0.0;
      double sum_alpha2_ylog = 0.0;
      for (int idx = 0; idx < mPi; ++idx) {
        double a2 = alpha2[idx];
        sum_lg_alpha2   += lgamma_cpp(a2);
        sum_alpha2_ylog += (a2 - 1.0) * y_log[idx];
      }
      loglik += lgamma_cpp(phi) - sum_lg_alpha2 + sum_alpha2_ylog;
      
      double S_phi_i = digamma_cpp(phi);
      double H_phi_i = psi1_phi;
      for (int idx = 0; idx < mPi; ++idx) {
        double mu2_i      = mu2[idx];
        double dig_alpha2 = digamma_cpp(alpha2[idx]);
        double tri_alpha2 = trigamma_cpp(alpha2[idx]);
        S_phi_i -= mu2_i * dig_alpha2;
        S_phi_i += mu2_i * y_log[idx];
        H_phi_i -= mu2_i * mu2_i * tri_alpha2;
      }
      S_phi += S_phi_i;
      H_phi += H_phi_i;
      
      // J: D x p
      arma::mat J(D, p, fill::zeros);
      arma::vec mu_sub(p);
      for (int k = 0; k < p; ++k)
        mu_sub[k] = mu_i[k + 1];
      
      for (int d = 0; d < D; ++d) {
        double mu_d = mu_i[d];
        for (int k = 0; k < p; ++k) {
          double indicator = (d == (k + 1)) ? 1.0 : 0.0;
          J(d, k) = mu_d * (indicator - mu_sub[k]);
        }
      }
      
      // Jp: mPi x p
      arma::mat Jp(mPi, p);
      for (int idx = 0; idx < mPi; ++idx)
        Jp.row(idx) = J.row(Pi[idx]);
      
      // Jsum: colsums(Jp)
      arma::rowvec Jsum = arma::sum(Jp, 0);
      
      // J2 = (Jp - outer(mu2, Jsum)) / Mi
      arma::mat J2(mPi, p);
      for (int r = 0; r < mPi; ++r)
        for (int k = 0; k < p; ++k)
          J2(r, k) = (Jp(r, k) - mu2[r] * Jsum[k]) / Mi;
      
      // h_i, psi1_alpha2, C
      arma::vec h_i(mPi);
      for (int idx = 0; idx < mPi; ++idx)
        h_i[idx] = y_log[idx] - digamma_cpp(alpha2[idx]);
      
      arma::vec psi1_alpha2(mPi);
      for (int idx = 0; idx < mPi; ++idx)
        psi1_alpha2[idx] = trigamma_cpp(alpha2[idx]);
      
      arma::mat C(mPi, mPi);
      C.fill(-psi1_phi);
      for (int r = 0; r < mPi; ++r)
        C(r, r) += psi1_alpha2[r];
      
      // JtCJ = crossprod(J2, C) %*% J2
      arma::mat JtCJ = J2.t() * C * J2;  // p x p
      
      // Jt_h = crossprod(J2, h_i)
      arma::vec Jt_h = J2.t() * h_i;     // p
      
      // S_vec update: S_vec <- S_vec + as.vector(phi * tcrossprod(x_i, Jt_h))
      // tcrossprod(x_i, Jt_h) is K x p, element (r,c) = x_i[r] * Jt_h[c]
      for (int k = 0; k < K; ++k) {
        double xik = x_i[k];
        for (int j = 0; j < p; ++j) {
          int idx = k + K * j; // column-major flattening
          S_vec[idx] += phi * xik * Jt_h[j];
        }
      }
      
      // I_mat update: I_mat <- I_mat + phi^2 * kronecker(JtCJ, txi[[i]])
      arma::mat txi_i = as<arma::mat>(txi[i]); // K x K
      for (int b = 0; b < p; ++b)
        for (int a = 0; a < p; ++a) {
          double coeff = phi * phi * JtCJ(a, b);
          for (int r = 0; r < K; ++r)
            for (int c = 0; c < K; ++c) {
              int row = r + K * a;
              int col = c + K * b;
              I_mat(row, col) += coeff * txi_i(r, c);
            }
        }
    } // end for i
    
    // penalized loglik
    double sum_beta2 = arma::accu(beta % beta);
    loglik_pen = loglik - 0.5 * lambda * sum_beta2;
    
    // beta_vec, S_vec_pen, I_mat_pen
    arma::vec beta_vec(dim);
    for (int j = 0; j < p; ++j)
      for (int k = 0; k < K; ++k)
        beta_vec[k + K * j] = beta(k, j);
    
    arma::vec S_vec_pen = S_vec - lambda_vec % beta_vec;
    arma::mat I_mat_pen = I_mat + lamI;
    
    if (std::abs(loglik_pen - loglik_old) < tol)
      break;
    
    // Newton step for beta: beta_vec_new <- beta_vec + solve(I_mat_pen, S_vec_pen)
    arma::vec delta = arma::solve(I_mat_pen, S_vec_pen);
    arma::vec beta_vec_new = beta_vec + delta;
    
    // reshape beta_vec_new back to beta (K x p, column-major)
    for (int j = 0; j < p; ++j)
      for (int k = 0; k < K; ++k)
        beta(k, j) = beta_vec_new[k + K * j];
    
    // Newton step for phi
    double phi_new = phi - S_phi / H_phi;
    if (phi_new <= 0.0) phi_new = 1e-4;
    phi = phi_new;
    
    loglik_old = loglik_pen;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("be")         = beta,
    Rcpp::Named("phi")        = phi,
    Rcpp::Named("loglik_pen") = loglik_pen
  );
}
