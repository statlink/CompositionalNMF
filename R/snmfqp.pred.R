# Define helper OUTSIDE the main function to avoid closure capture
.solve_w_row_qp <- function(i, G_W, g_W, A_W, b_W) {
  sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
  abs(sol$solution)
}

snmfqp.pred <- function(xnew, H, ridge = 1e-8, ncores = 1) {

  runtime <- proc.time()
  n <- dim(xnew)[1]  ;  k <- dim(H)[1]

  A_W <- t( rbind(rep(1, k), diag(k)) )
  b_W <- c(1, rep(0, k))
  ridgek <- diag(ridge, k)
  G_W <- 2 * tcrossprod(H) + ridgek
  g_W <- 2 * tcrossprod(H, xnew)

  suppressWarnings({

  # CREATE CLUSTER ONCE BEFORE LOOP (like snmf.sqp)
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(quadprog))
    parallel::clusterExport( cl, varlist = c(".solve_w_row_qp", "A_W", "b_W"), envir = environment() )
  }

  if ( ncores > 1 ) {
    # Export iteration-specific variables
    parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment())
      Wnew <- t( parallel::parSapply(cl, 1:n, .solve_w_row_qp, G_W = G_W, g_W = g_W, A_W = A_W, b_W = b_W) )
  } else {
    Wnew <- matrix(nrow = n, ncol = k)
    for ( i in 1:n ) {
      sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
      Wnew[i, ] <- abs(sol$solution)
    }
  }
  Znew <- Wnew %*% H

  })  # end suppressWarnings

  runtime <- proc.time() - runtime

  list(Wnew = Wnew, Znew = Znew, runtime = runtime)
}


