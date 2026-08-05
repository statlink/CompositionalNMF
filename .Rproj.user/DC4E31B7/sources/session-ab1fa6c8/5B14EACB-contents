\name{snmfqp.pred}
\alias{snmfqp.pred}
\title{
Prediction of W and Z for new simplicial data
}

\description{
Prediction of \eqn{W} and \eqn{Z} matrices for new simplicial data.
}

\usage{
snmfqp.pred(xnew, H, ridge = 1e-8, ncores = 1)
}

\arguments{
\item{xnew}{
An \eqn{n \times D} matrix with some new compositional data. Zero values are allowed.
}
\item{H}{
The estimated H from the \code{link{snmf.qp}} function.
}
\item{ridge}{
A small quantity added in the diagonal of the \eqn{D} matrix.
}
\item{ncores}{
Do you want the estimation to be performed in parallel? If yes, specify the number of cores to use.
}
}

\details{
Based on an already SNMF that was produced by minimizing the square of the Frobenius norm, the function
estimates the \eqn{W} and \eqn{Z} matrices for some new data.
}

\value{
\item{Wnew}{
The \eqn{W} matrix, an \eqn{n \times k} matrix with the mapped compositional data.
}
\item{Znew}{
The reconstructed data, \eqn{Znew = WnewH}.
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
mod <- snmf.qp(x, 3)
est <- snmfqp.pred(x[1:10, ], mod$H)
}

