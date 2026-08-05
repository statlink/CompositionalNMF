\name{snmf.sqp}
\alias{snmf.sqp}
\title{
Simplicial NMF minimizing the Frobenius norm with sparse data
}

\description{
Simplicial NMF minimizing the Frobenius norm with sparse data.
}

\usage{
snmf.sqp(x, k, W = NULL, H = NULL, maxiter = 1000, tol = 1e-4, ridge = 1e-8,
history = FALSE, ncores = 1)
}

\arguments{
\item{x}{
An \eqn{n \times D} dgC class sparse matrix with compositional data.
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
The tolerance value to terminate the quadratic programming algorithm. The value is set to 1e-4 in this case because with large scale and/or sparse data the computation time is really high. So, we sacrifice some accuracy over speed.
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
Simplicial nonnegative matrix factorization using quadratic programming is performed. This is the classical NMF with simplicial constraints on the matrices \eqn{W} and \eqn{H}. The objective function to be minimized is the square of the Frobenius norm. This function is suitable for large scale sparse data, and parallel computing is a must in this case. Note that we do not use k-means here and that the reconstruced matrix Z is not returned with this function for capacity purposes.
}

\value{
\item{W}{
The \eqn{W} matrix, an \eqn{n \times k} matrix with the mapped compositional data.
}
\item{H}{
The \eqn{H} matrix, an \eqn{k \times D} matrix. Each row of this matrix is a compositional vector.
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
Alenazi A. and Tsagris M. (2026). Simplicial nonnegative matrix factorization. In preparation.

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
group <- as.numeric(iris[, 5])
ternary(mod$W, colour = group)
}

