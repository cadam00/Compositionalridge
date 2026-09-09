\name{Ridge alpha-regression for compositional data}
\alias{areg.ridge}
\title{
Ridge \eqn{\alpha}-regression for compositional data
}

\description{
Ridge \eqn{\alpha}-regression for compositional data
}

\usage{
areg.ridge(y, x, a, lambda = NULL, nlambda = 100, yb = NULL, xnew = NULL)
}

\arguments{
\item{y}{
A matrix with the compositional data.
}
\item{x}{
A matrix with the continuous predictor variables or a data frame including categorical predictor variables.
}
\item{a}{
The values of \eqn{\alpha} for the \eqn{\alpha}-transformation. If zero values are present
only the positive values of \eqn{\alpha} are used.
}
\item{lambda}{
A vector with a sequence of the ridge \eqn{\lambda} values. If you do not know it leave it NULL
and the function will compute it internally.
}
\item{nlambda}{
The number of \eqn{\lambda} values to compute.
}
\item{yb}{If you have already \eqn{\alpha}-transformed the data using \eqn{\alpha} as given
in the argument "a", put it here. Othewrise leave it NULL.
}
\item{xnew}{
If you have new data whose response values you want to compute, pass them here.
}
}

\details{
The ridge $\eqn{\alpha}-regression is applied.
}

\value{
A list including:
\item{runtime}{
The runtime of the ridge \eqn{\alpha}-regression.
}
\item{info}{
A matrix with two columns containing the lambda, and the Kullback-Leibler
divergence between the observed and the fitted compositions.
}
\item{be}{
A list with the ridge beta coefficients.
}
\item{est}{
A list with the ridge estimated compositions, if xnew is not NULL.
}
}

\references{
Tsagris M. and Pantazis Y. (2026). The \eqn{\alpha}--regression for compositional data:
a unified framework for standard, temporal and spatial regression models including compositional predictors.
https://arxiv.org/pdf/2510.12663

Tsagris M. (2015). Regression analysis with compositional data containing zero values.
Chilean Journal of Statistics, 6(2): 47-57.
https://arxiv.org/pdf/1508.01913v1.pdf

Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for
compositional data.
In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{cv.aregridge}, \link{zadr.ridge}
}
}

\examples{
x <- matrix( rnorm(150 * 10), ncol = 10 )
y <- as.matrix(iris[, 1:3])
y <- y / rowSums(y)
ind <- sample(150, 10)
for ( k in ind )  y[k, sample(3, 1)] <- 0
y <- y / rowSums(y)
mod <- areg.ridge(y, x, 0.5)
}

