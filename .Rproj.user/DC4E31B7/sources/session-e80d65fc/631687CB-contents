\name{snmf.eg}
\alias{snmf.eg}
\title{
Simplicial NMF minimizing the Frobenius norm
}

\description{
Simplicial NMF minimizing the Frobenius norm using exponentiated gradient descent.
}

\usage{
snmf.eg(x, k, W = NULL, H = NULL, lr_w = 0.1, lr_h = 0.1, maxiter = 10000, tol = 1e-6)
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
\item{lr_w}{
The learning rate to multiply the gradient for the \eqn{W} matrix.
}
\item{lr_h}{
The learning rate to multiply the gradient for the \eqn{H} matrix.
}
\item{maxiter}{
The maximum number of iterations allowed.
}
\item{tol}{
The tolerance value to terminate the exponentiated gradient descent algorithm.
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
\code{link{ternary}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- snmf.eg(x, 3)
group <- as.numeric(iris[, 5])
ternary(mod$W, colour = group)
}

