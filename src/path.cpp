#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
Rcpp::List path(
    arma::mat beta,
    double phi,
    const arma::mat& ly,
    const arma::mat& x,
    double lambda,
    const arma::vec& pen,
    int maxit,
    double tol
) {
  const int n = ly.n_rows;
  const int D = ly.n_cols;
  const int p = D - 1;
  const int K = x.n_cols;
  const int dim = K * p;
  
  // Match the original .path() initialization.
  phi = 1.0;
  
  // Ridge penalty: pen[0] = 0 leaves the intercept unpenalized.
  arma::vec lambda_vec(dim);
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < K; ++k) {
      lambda_vec[k + K * j] = lambda * pen[k];
    }
  }
  
  arma::mat lamI = arma::diagmat(lambda_vec);
  
  // Column-major vectorization: beta(k,j) -> k + K*j.
  arma::vec beta_vec(dim);
  for (int j = 0; j < p; ++j) {
    for (int k = 0; k < K; ++k) {
      beta_vec[k + K * j] = beta(k, j);
    }
  }
  
  double loglik_old = -std::numeric_limits<double>::infinity();
  double loglik_pen = NA_REAL;
  
  arma::mat eta(n, p);
  arma::mat exp_eta(n, p);
  arma::mat mu(n, D);
  arma::mat alpha(n, D);
  
  for (int iter = 0; iter < maxit; ++iter) {
    
    // Multinomial mean.
    eta = x * beta;
    exp_eta = arma::exp(eta);
    
    arma::vec sum_exp = arma::sum(exp_eta, 1);
    
    mu.col(0) = 1.0 / (1.0 + sum_exp);
    
    for (int j = 0; j < p; ++j) {
      mu.col(j + 1) =
        exp_eta.col(j) / (1.0 + sum_exp);
    }
    
    alpha = phi * mu;
    
    const double psi1_phi = R::trigamma(phi);
    
    double loglik = 0.0;
    double S_phi = 0.0;
    double H_phi = 0.0;
    
    arma::vec S_vec(dim, fill::zeros);
    arma::mat I_mat(dim, dim, fill::zeros);
    
    // Accumulate observation-level contributions.
    for (int i = 0; i < n; ++i) {
      
      const arma::rowvec x_i = x.row(i);
      const arma::rowvec mu_i = mu.row(i);
      const arma::rowvec alpha_i = alpha.row(i);
      const arma::rowvec ly_i = ly.row(i);
      
      arma::rowvec h(D);
      arma::rowvec tri_alpha(D);
      
      double loglik_i = R::lgammafn(phi);
      double S_phi_i = R::digamma(phi);
      double H_phi_i = psi1_phi;
      
      for (int d = 0; d < D; ++d) {
        const double a = alpha_i[d];
        const double m = mu_i[d];
        const double l = ly_i[d];
        
        loglik_i -= R::lgammafn(a);
        loglik_i += (a - 1.0) * l;
        
        h[d] = l - R::digamma(a);
        tri_alpha[d] = R::trigamma(a);
        
        S_phi_i += m * h[d];
        H_phi_i -= m * m * tri_alpha[d];
      }
      
      loglik += loglik_i;
      S_phi += S_phi_i;
      H_phi += H_phi_i;
      
      // Jacobian: J(d,j) = d mu_d / d eta_j.
      arma::mat J(D, p);
      
      for (int d = 0; d < D; ++d) {
        for (int j = 0; j < p; ++j) {
          const double indicator = (d == j + 1) ? 1.0 : 0.0;
          
          J(d, j) =
            mu_i[d] * (indicator - mu_i[j + 1]);
        }
      }
      
      // Calculate J' h and J' C J.
      arma::vec Jt_h(p, fill::zeros);
      arma::mat JtCJ(p, p, fill::zeros);
      arma::rowvec Jsum = arma::sum(J, 0);
      
      for (int d = 0; d < D; ++d) {
        for (int a = 0; a < p; ++a) {
          Jt_h[a] += J(d, a) * h[d];
          
          for (int b = 0; b < p; ++b) {
            JtCJ(a, b) +=
              J(d, a) * tri_alpha[d] * J(d, b);
          }
        }
      }
      
      // C = diag(trigamma(alpha)) -
      //     trigamma(phi) * 1 * 1'.
      JtCJ -= psi1_phi * Jsum.t() * Jsum;
      
      // Coefficient score.
      for (int j = 0; j < p; ++j) {
        for (int k = 0; k < K; ++k) {
          S_vec[k + K * j] +=
            phi * x_i[k] * Jt_h[j];
        }
      }
      
      // Information matrix, equivalent to
      // phi^2 * kronecker(JtCJ, tcrossprod(x_i)).
      for (int a = 0; a < p; ++a) {
        for (int b = 0; b < p; ++b) {
          const double coeff = phi * phi * JtCJ(a, b);
          
          for (int r = 0; r < K; ++r) {
            for (int c = 0; c < K; ++c) {
              I_mat(r + K * a, c + K * b) +=
                coeff * x_i[r] * x_i[c];
            }
          }
        }
      }
    }
    
    // Correct penalized log-likelihood.
    // Uses lambda_vec, so the intercept remains unpenalized
    // whenever pen[0] = 0.
    loglik_pen =
      loglik - 0.5 * arma::dot(lambda_vec, beta_vec % beta_vec);
    
    // Penalized score and information.
    arma::vec S_vec_pen = S_vec - lambda_vec % beta_vec;
    arma::mat I_mat_pen = I_mat + lamI;
    
    // Preserve the original convergence check.
    if (std::abs(loglik_pen - loglik_old) < tol) {
      break;
    }
    
    // Newton update for beta.
    arma::vec delta = arma::solve(I_mat_pen, S_vec_pen);
    beta_vec += delta;
    
    for (int j = 0; j < p; ++j) {
      for (int k = 0; k < K; ++k) {
        beta(k, j) = beta_vec[k + K * j];
      }
    }
    
    // Newton update for phi, matching the original R behavior.
    double phi_new = phi - S_phi / H_phi;
    
    if (phi_new <= 0.0) {
      phi_new = 1e-4;
    }
    
    phi = phi_new;
    loglik_old = loglik_pen;
  }
  
  return Rcpp::List::create(
    Rcpp::Named("be") = beta,
    Rcpp::Named("phi") = phi,
    Rcpp::Named("loglik_pen") = loglik_pen
  );
}