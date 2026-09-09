\name{Cross-validation for ridge regression for compositional data}
\alias{cv.zadrridge}
\alias{cv.drridge}
\alias{cv.kldridge}
\alias{cv.mzilnridge}
\title{
Cross-validation for ridge regression for compositional data
}

\description{
Cross-validation for ridge regression for compositional data.
}

\usage{
cv.zadrridge(y, x, lambda = NULL, nlambda = 100, tol = 1e-06, maxit = 100,
folds = NULL, nfolds = 10, seed = NULL)
cv.drridge(y, x, lambda = NULL, nlambda = 100, tol = 1e-06, maxit = 100,
folds = NULL, nfolds = 10, seed = NULL)
cv.kldridge(y, x, lambda = NULL, nlambda = 100,
folds = NULL,  nfolds = 10, seed = NULL)
cv.mzilnridge(y, x, lambda = NULL, nlambda = 100,
folds = NULL, nfolds = 10, seed = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data (dependent variable).
Only for the Dirichlet regression, zeros are not allowed.
}
\item{x}{
The predictor variable(s), they can only be continnuous variables.
}
\item{lambda}{
A vector with a sequence of the ridge \eqn{\lambda} values. If you do not know it leave it NULL
and the function will compute it internally.
}
\item{nlambda}{
The number of \eqn{\lambda} values to compute.
}
\item{tol}{
The tolerance value to terminate the Fisher scoring algorithm.
}
\item{maxit}{
The maximum number of iterations the Fisher scoring algorithm will perform.
}
\item{folds}{
If you have the list with the folds supply it here. You can also leave it NULL and it will create folds.
}
\item{nfolds}{
The number of folds to produce.
}
\item{seed}{
You can specify your own seed number here or leave it NULL.
}
}

\details{
K-fold cross-validation for the ridge regression using the Dirichlet regression, the
zero adjuster Dirichlet regression, the Kullback-Leilber divergence regression or
the multivariate zero inflated logistic normal regression.
}

\value{
A list including:
\item{runtime}{
The runtime of the ridge regression model.
}
\item{info}{
A matrix with two or three columns containing the lambda, the penalized log-likelihood (for dr and zadr),
and the Kullback-Leibler divergence between the observed and te fitted compositions.
}
\item{be}{
A list with the ridge beta coefficients.
}
\item{sest}{
A list with the ridge estimated compositions.
}
}

\references{
Li Z., Lee K., Karagas M. R., Madan J. C., Hoen A. G., O'Malley A. J. and Li H. (2018).
Conditional regression based on a multivariate zero-inflated logistic-normal model for
microbiome relative abundance data.
Statistics in Biosciences, 10(3): 587--608.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.zadrridge}, \link{areg.ridge}
}
}

\examples{
x <- matrix( rnorm(150 * 10), ncol = 10 )
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
mod <- kld.ridge(y, x)
}

