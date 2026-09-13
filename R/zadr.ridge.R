zadr.ridge <- function(y, x, lambda = NULL, nlambda = 100, xnew = NULL, tol = 1e-6, maxit = 100) {

  runtime <- proc.time()
  n <- dim(y)[1]   ;  D <- dim(y)[2]  ;  p <- D - 1
  x <- as.matrix(x)
  K <- dim(x)[2] + 1
  nam <- colnames(x)
  if ( !is.null(nam) ) {
    nam <- c( "Intercept", nam )
  } else   nam <- c( "Intercept", paste("X", 1:(K - 1), sep = "" ) )
  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x, std = TRUE)
  x <- t( ( t(x) - m ) / s )
  x <- cbind(1, x)
  logy <- log(y)
  txi <- list()
  for ( i in 1:n )  txi[[ i ]] <- tcrossprod(x[i, ])
  pen <- rep(1, K)
  pen[1] <- 0
  mod <- kld.ridge(y, x[, -1], lambda = lambda, nlambda = nlambda)
  lambda <- mod$info[, 1]
  nlambda <- length(lambda)
  ini <- mod$be

  pos <- y > 0
  a1 <- which( Rfast::rowsums( pos ) == D )
  a2 <- which( Rfast::rowsums( pos ) != D )
  n1 <- length(a1)   ;   n2 <- n - n1
  if ( n2 > 0 ) {
    za <- y[a2, , drop = FALSE]
    za[za == 0] <- 1  ;  za[za < 1] <- 0
    theta <- as.vector( table( apply(za, 1, paste, collapse = ",") ) )
    const <- n1 * log(n1 / n) + sum( theta * log(theta / n) )
  } else  const <- n1 * log(n1 / n)
  Pidx <- lapply( 1:n, function(i)  which( pos[i, ] ) )

  be <- list()   ;   loglik_pen <- phi <- numeric(nlambda)
  beta <- ini[[ 1 ]][, -1] - ini[[ 1 ]][, 1]

  for ( vim in 1:nlambda ) {
    mod <- zpath(beta, 1.0, logy, y, x, txi, Pidx, lambda[vim], n, D, p, K, pen, maxit, tol)
    beta <- mod$be
    colnames(beta) <- paste0("Y", 2:D)
    rownames(beta) <- nam
    be[[ vim ]] <- beta
    phi[vim] <-  mod$phi
    loglik_pen[vim] <- mod$loglik_pen + const
  }

  est <- NULL
  if ( !is.null(xnew) ) {
    est <- list()
    xnew <- as.matrix(xnew)
    if ( nrow(xnew) == 1 )  xnew <- t(xnew)
    xnew <- t( ( t(xnew) - m ) / s )
    xnew <- cbind(1, xnew)
    for ( i in 1:nlambda ) {
      mu <- xnew %*% be[[ i ]]
      mu <- cbind(1, exp(mu))
      est[[ i ]] <- mu / Rfast::rowsums(mu)
    }
    names(est) <- paste("lambda=", lambda, sep="")
  }
  mu <- NULL

  info <- cbind(lambda, loglik_pen, phi)
  colnames(info) <- c("lambda", "loglik", "phi")
  runtime <- proc.time() - runtime

  list( runtime = runtime, info = info, be = be, est = est )
}




# .zpath <- function(beta, phi, logy, y, x, txi, Pidx, lambda, n, D, p, K, pen, maxit, tol) {
# 
#   lambda_vec <- lambda * rep(pen, p)
#   lamI <- diag( lambda_vec, nrow = length(lambda_vec) )
#   phi <- 1.0
#   loglik_old <-  -Inf
# 
#   for ( iter in 1:maxit ) {
#     eta <- x %*% beta
#     exp_eta <- exp(eta)
#     mu <- cbind(1, exp_eta) / ( 1 + Rfast::rowsums(exp_eta) )
#     psi1_phi <- trigamma(phi)
#     loglik <- S_vec <- I_mat <- S_phi <- H_phi <- 0
# 
#     for ( i in 1:n ) {
#       Pi <- Pidx[[ i ]]
#       x_i <- x[i, ]
#       mu_i <- mu[i, ]
#       Mi <- sum( mu_i[Pi] )
#       mu2 <- mu_i[Pi] / Mi
#       alpha2 <- phi * mu2
#       y_log <- logy[i, Pi]
# 
#       loglik <- loglik + lgamma(phi) - sum( lgamma(alpha2) ) + sum( (alpha2 - 1) * y_log )
#       S_phi <- S_phi + digamma(phi) - sum( mu2 * digamma(alpha2) ) + sum( mu2 * y_log )
#       H_phi <- H_phi + trigamma(phi) - sum( mu2^2 * trigamma(alpha2) )
#       J <- matrix(0, nrow = D, ncol = p)
#       mu_sub <- mu_i[2:D]
#       for ( d in 1:D ) {
#         for ( k in 1:p ) {
#           J[d, k] <- mu_i[d] * ( ( d == (k + 1) ) - mu_sub[k] )
#         }
#       }
#       Jp <- J[Pi, , drop = FALSE]
#       Jsum <- Rfast::colsums(Jp)
#       J2 <- ( Jp - outer(mu2, Jsum) ) / Mi
#       h_i <- y_log - digamma(alpha2)
#       psi1_alpha2 <- trigamma(alpha2)
#       C <- diag( psi1_alpha2, nrow = length(Pi) ) - psi1_phi
#       JtCJ <- crossprod(J2, C) %*% J2
#       Jt_h <- crossprod(J2, h_i)
#       S_vec <- S_vec + as.vector( phi * tcrossprod(x_i, Jt_h) )
#       I_mat <- I_mat + phi^2 * kronecker(JtCJ, txi[[ i ]] )
#     }  ##  end for (i in 1:n)
# 
#     loglik_pen <- loglik - 0.5 * lambda * sum( beta^2 )
#     beta_vec <- as.vector(beta)
#     S_vec_pen <- S_vec - lambda_vec * beta_vec
#     I_mat_pen <- I_mat + lamI
# 
#     if ( abs(loglik_pen - loglik_old) < tol ) break
# 
#     beta_vec_new <- beta_vec + solve(I_mat_pen, S_vec_pen)
#     beta <- matrix(beta_vec_new, nrow = K, ncol = p)
#     phi_new <- phi - S_phi / H_phi
#     if ( phi_new <= 0 )  phi_new <- 1e-4
#     phi <- phi_new
#     loglik_old <- loglik_pen
#   }
#   list(be = beta, phi = phi, loglik_pen = loglik_pen)
# }
