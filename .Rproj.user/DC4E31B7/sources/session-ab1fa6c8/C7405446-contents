\name{snmfqp.cv}
\alias{snmfqp.cv}
\title{
K-fold cross-validation for choosing the rank in SNMF
}

\description{
K-fold cross-validation for choosing the rank in SNMF.
}

\usage{
snmfqp.cv(x, k = 3:10, k_means = TRUE, veo = FALSE, lr_h = 0.1, maxiter = 1000,
tol = 1e-6, ridge = 1e-8, ncores = 1, folds = NULL, nfolds = 10, graph = FALSE)
}

\arguments{
\item{x}{
An \eqn{n \times D} matrix with compositional data. Zero values are allowed.
}
\item{k}{
The number of lower dimensions. It must be less than the dimensionality of the data, at most \eqn{D-1}.
}
\item{k_means}{
If this is TRUE, then the K-means algorithm is used to initiate the W and H matrices.
}
\item{veo}{
If the number of variables excceeds the number of observations set this is equal to TRUE.
In this case, the sparse k-means algorithm of Witten and Tibshirani (2010) is used to initialize the H matrix.
}
\item{lr_h}{
If veo is TRUE, then the exponentiated gradient descent method is used to update the H matrix. In this case you need to supply the value of the learning rate, which is 0.1 by default.
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
\item{folds}{
If you have the list with the folds supply it here.
You can also leave it NULL and it will create folds.
}
\item{nfolds}{
The number of folds to produce.
}
\item{graph}{
If this is TRUE, the plot of the predicted error will be plotted.
}
}

\details{
K-fold cross-validation to select the optimal rank k.
}

\value{
\item{sse}{
The matrix with the sum of squares of residuals.
}
\item{mspe}{
A vector with the mean squares of residuals.
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
\code{ \link{snmf.qp}, \link{snmfqp.pred}, link{ternary}, \link{pyramid}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
x <- x / rowSums(x)
mod <- snmf.qp(x, 3)
est <- snmfqp.pred(x[1:10, ], mod$H)
}

