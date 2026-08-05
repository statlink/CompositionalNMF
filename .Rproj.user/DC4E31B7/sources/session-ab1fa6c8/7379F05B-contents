snmfqp.cv <- function(x, k = 3:10, k_means = TRUE, veo = FALSE, lr_h = 0.1, maxiter = 1000,
             tol = 1e-6, ridge = 1e-8, ncores = 1, folds = NULL, nfolds = 10, graph = FALSE) {

  runtime <- proc.time()
  n <- dim(x)[1]
  if ( is.null(folds) )  folds <- Compositional::makefolds(1:n, nfolds = 10, stratified = FALSE)
  K <- length(folds)
  sse <- matrix( nrow = K, ncol = length(k) )
  rownames(sse) <- paste("Fold ", 1:K, sep = "")
  colnames(sse) <- paste("k=", k, sep = "")

  for ( i in 1:K ) {
    xtrain <- x[ -folds[[ i ]], ]
    xtest <- x[ folds[[ i ]], ]
    for ( j in 1:length(k) ) {
      mod <- CompositionalNMF::snmf.qp(xtrain, k[j])
      test <- CompositionalNMF::snmfqp.pred(xtest, mod$H)$Znew
      sse[i, j] <- sum( ( xtest - test)^2 )
    }
  }

  mspe <- Rfast::colmeans(sse)
  if ( graph ) {
    plot(mspe, type = "b", pch = 16, col = 4, lwd = 2, xlab = "Rank (k)",
         ylab = "Frobenius norm", xaxt = "n", cex.lab = 1.3, cex.axis = 1.3)
    axis(1, at = 1:length(k), labels = k)
  }

  runtime <- proc.time() - runtime
  list( sse = sse, mspe = mspe, runtime = runtime )
}
