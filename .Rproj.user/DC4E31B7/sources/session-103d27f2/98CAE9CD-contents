areg.ridge <- function(y, x, a, lambda = NULL, nlambda = 100, yb = NULL, xnew = NULL) {

  reg <- function(para, ya, ax, a, ha, d, D, lambda, pen) {
    be <- matrix(para, ncol = d)
    zz <- cbind( 1, exp(ax %*% be) )
    ta <- rowSums(zz)
    za <- zz / ta
    ma <- ( D / a * za - 1/a ) %*% ha
    r <- as.vector(ya - ma)
    rp <- sqrt(0.5 * lambda) * para[pen]
    c(r, rp)
  }

  ## Analytical Jacobian of reg()
  jac <- function(para, ya, ax, a, ha, d, D, lambda, pen) {
    n <- dim(ax)[1]   ;   p <- dim(ax)[2]
    beta_idx  <- lapply(1:d, function(k) ( (k - 1L) * p + 1L):(k * p) )
    resid_idx <- lapply(1:d, function(j) ( (j - 1L) * n + 1L):(j * n) )
    be <- matrix(para, ncol = d)
    zz <- cbind(1, exp(ax %*% be))
    za <- zz / rowSums(zz)
    zh <- za %*% ha
    J  <- matrix(0, nrow = n * d, ncol = d * p)
    Da <- D / a
    for ( j in 1:d ) {
      for ( k in 1:d ) {
        w <- Da * za[, k + 1L] * (ha[k + 1L, j] - zh[, j])
        J[ resid_idx[[j]], beta_idx[[k]] ] <-  - w * ax
      }
    }
    npen <- length(pen)
    Jp <- matrix(0, nrow = npen, ncol = length(para))
    Jp[ cbind(seq_len(npen), pen) ] <- sqrt(0.5 * lambda)
    rbind(J, Jp)
  }

  runtime <- proc.time()

  m <- Rfast::colmeans(x)
  s <- Rfast::colVars(x)
  x <- t( ( t(x) - m ) / s )

  if ( is.null(yb) ) {
    ya <- Compositional::alfa(y, a)$aff
  } else  ya <- yb
  x <- cbind(1, x)
  ax <- a * x

  D <- dim(y)[2]
  d <- D - 1  ## dimensionality of the simplex
  p <- dim(ax)[2]
  intercept_idx <- ( 1:d - 1 ) * p + 1
  pen <- setdiff( seq_len(d * p), intercept_idx)
  npen <- length(pen)

  if ( is.null( colnames(x) ) ) {
    nama <- c("constant", paste("X", 1:(p - 1), sep = "") )
  } else  nama <- c("constant", colnames(x)[-1] )

  if ( a <= 1e-5 ) {
    mod <- Compositional::comp.reg(y, x[, -1], yb = yb)
    be <- mod$be

  } else {
    ha  <- t( Compositional::helm(D) )
    ridge <- glmnet::glmnet(x[, -1], ya, alpha = 0, lambda = lambda, nlambda = nlambda,
                            standardize = FALSE, family = "mgaussian")
    lambda <- ridge$lambda
    nlambda <- length(ridge$lambda)

    be <- list()
    dev <- numeric(nlambda)
    for ( i in 1:nlambda ) {
      ini <- NULL
      for ( j in 1:d )  ini <- c(ini, coef(ridge)[[ j ]][, i])
        suppressWarnings({
        mod <- minpack.lm::nls.lm(par = ini, fn = reg, jac = jac, ya = ya, ax = ax, a = a,
               ha = ha, d = d, D = D, lambda = lambda[i], pen = pen)
      })
      be[[ i ]] <- matrix(mod$par, ncol = d)
      rownames(be[[ i ]]) <- nama
      colnames(be[[ i ]]) <- paste0("Y", 2:D)
      dev[i] <- mod$deviance
    }
  }  ## end if (a == 0)

  names(be) <- paste("lambda=", lambda, sep="")
  info <- cbind(lambda, dev)
  colnames(info) <- c("lambda", "deviance")

  # prediction
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
  runtime <- proc.time() - runtime

  list(runtime = runtime, be = be, info = info, est = est)
}
