\name{snmf.qp}
\alias{snmf.qp}
\title{
Simplicial NMF minimizing the Frobenius norm
}

\description{
Simplicial NMF minimizing the Frobenius norm using quadratic programming.
}

\usage{
snmf.qp(x, k, W = NULL, H = NULL, k_means = TRUE, lr_h = 0.1,
maxiter = 1000, tol = 1e-6, ridge = 1e-8, history = FALSE, ncores = 1)
}

\arguments{
\item{x}{
An \eqn{n \times D} matrix with compositional data. Zero values are allowed.
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
\item{k_means}{
If this is TRUE, then the K-means algorithm is used to initiate the W and H matrices.
}
\item{lr_h}{
If the number of variables excceeds the number of observations, the sparse k-means algorithm of
Witten and Tibshirani (2010) is used to initialize the H matrix. The exponentiated gradient descent
method is used to update the H matrix. In this case you need to supply the value of the learning rate,
which is 0.1 by default.
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
\item{history}{
If this is TRUE, the reconstruction error at each iteration is returned.
}
\item{ncores}{
Do you want the update of W to be performed in parallel? If yes, specify the number of cores to use.
}
}

\details{
Simplicial nonnegative matrix factorization using quadratic programming is performed. This is the classical NMF with simplicial constraints on the matrices \eqn{W} and \eqn{H}. The objective function to be minimized is the square of the Frobenius norm.
}

\value{
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
\item{error}{
If the argument history was set to TRUE the reconstruction error at each iteration will be performed, otherwise this is NULL.
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
\code{ \link{snmf.sqp}, link{ternary}, \link{pyramid}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- snmf.qp(x, 3)
group <- as.numeric(iris[, 5])
ternary(mod$W, colour = group)
}

