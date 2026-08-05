\name{snmfqp.reg}
\alias{snmfqp.reg}
\title{
Simplicial NMF minimizing the Frobenius norm
}

\description{
Simplicial NMF minimizing the Frobenius norm using quadratic programming.
}

\usage{
snmfqp.reg(x, z, k, W = NULL, H = NULL, maxiter = 1000, tol = 1e-6,
ridge = 1e-8, ncores = 1)
}

\arguments{
\item{x}{
An \eqn{n \times D} matrix with compositional data. Zero values are allowed.
}
\item{z}{
An \eqn{n \times D} matrix with the compositional predictors. Zero values are allowed.
}
\item{k}{
The number of lower dimensions. It must be less than the dimensionality of the data, at most \eqn{D-1}.
}
\item{W}{
If you have an initial estimate for W supply it here. Otherwise leave it NULL.
}
\item{H}{
If you have an initial estimate for H supply it here, otherwise leave it NULL.
}
\item{maxiter}{
The maximum number of iterations allowed.
}
\item{tol}{
The tolerance value to terminate the quadratic programming algorithm.
}
\item{ridge}{
A small quantity added in the diagonal of the \eqn{D} matrix.
}
\item{ncores}{
Do you want the update of W to be performed in parallel? If yes, specify the number of cores to use.
}
}

\details{
Simplicial nonnegative matrix factorization using quadratic programming is performed.
But now, a simplicial predictor is used. So this is SNMF regression style.
}

\value{
\item{B}{
The \eqn{B} compositional matrix of coefficients of \eqn{z}.
}
\item{W}{
The \eqn{W} matrix, an \eqn{n \times k} matrix with the mapped compositional data.
}
\item{H}{
The \eqn{H} matrix, an \eqn{k \times D} matrix. Each row of this matrix is a compositional vector.
}
\item{Z}{
The reconstructed data, \eqn{Z = WH}.
}
\item{obj}{
The reconstruction error, \eqn{ ||x - Z||_F^2}.
}
\item{iters}{
The number of iterations performed.
}
\item{runtime}{
The runtime required by the algorithm.
}
}

\references{
Alenazi A. and Tsagris M. (2026). Simplicial nonnegative matrix factorization with applications
to high--dimensional and large--scale datasets. In preparation.

Witten D. M. and Tibshirani R. (2010). A framework for feature selection in clustering.
Journal of the American Statistical Association, 105(490): 713--726.

Cutler A. and Breiman L. (1994). Archetypal analysis. Technometrics, 36(4): 338--347.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{snmf.qp}, link{ternary}, \link{pyramid}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
z <- matrix(runif(150 * 5), ncol = 5)
z <- z / rowSums(z)
mod <- snmfqp.reg(x, z, 3)
}

