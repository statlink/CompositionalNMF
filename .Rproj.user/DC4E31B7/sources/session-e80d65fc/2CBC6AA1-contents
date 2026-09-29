snmfqp.reg <- function(x, z, k, W = NULL, H = NULL, maxiter = 1000, tol = 1e-6, ridge = 1e-8, ncores = 1) {
  n <- dim(x)[1]
  D <- dim(x)[2]  # number of compositional parts
  q <- dim(z)[2]  # number of covariates

  runtime <- proc.time()

  if ( is.null(W) ) {
    # Initialize with simplex constraints
    W <- matrix( rangen::Runif(n * k), nrow = n, ncol = k )
    W <- W / Rfast::rowsums(W)
  }
  if ( is.null(H) ) {
    H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D )
    H <- H / Rfast::rowsums(H)
  }
  B <- matrix(nrow = q, ncol = D)

  # Setup QP constraints for simplex (rows sum to 1, non-negative)
  # For B: each row sums to 1
  A_B <- t( rbind(rep(1, D), diag(D)) )
  b_B <- c(1, rep(0, D))

  # For W: each row sums to 1
  A_W <- t( rbind(rep(1, k), diag(k)) )
  b_W <- c(1, rep(0, k))
  ridgek <- diag(ridge, k)
  # For H: each row sums to 1 (D-dimensional simplex)
  A_H <- A_B  # Same as B: D-dimensional simplex
  b_H <- b_B
  diagD <- diag(D)

  sse_old <- Inf
  cl <- NULL

  # Setup parallel cluster if needed
  if ( ncores > 1 ) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(quadprog))
  }

  for ( iter in 1:maxiter ) {
    ## Update B ( q x D, each row is a simplex )
    R <- x - W %*% H  # residual after removing W*H
    for ( i in 1:q ) {
      # For row i: minimize ||R - z[,i] %*% t(B[i,])||
      zi_norm_sq <- sum(z[, i]^2)
      G_B <- (2 * zi_norm_sq + ridge) * diagD
      dvec <- drop( 2 * crossprod(z[, i], R) )
      sol <- quadprog::solve.QP(Dmat = G_B, dvec = dvec, Amat = A_B, bvec = b_B, meq = 1)
      B[i, ] <- abs( sol$solution )
    }

    ## Update W ( n x k, each row is a simplex )
    R <- x - z %*% B  # residual after removing z*B
    G_W <- 2 * tcrossprod(H) + ridgek
    g_W <- 2 * tcrossprod(H, R)
    if ( ncores > 1 ) {
      parallel::clusterExport( cl, varlist = c("G_W", "g_W", "A_W", "b_W"),
                             envir = environment() )
      W <- t( parallel::parSapply(cl, 1:n, function(i) {
           sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
           abs( sol$solution )
      }) )
    } else {
      for ( i in 1:n ) {
        sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
        W[i, ] <- abs( sol$solution )
      }
    }

    ## Update H ( k x D, each row is a simplex )
    # R already computed above (x - z*B)
    for ( i in 1:k ) {
      # For row i: minimize ||R - W[,i] %*% t(H[i,])||
      wi_norm_sq <- sum(W[, i]^2)
      G_H <- ( 2 * wi_norm_sq + ridge ) * diagD
      dvec <- drop( 2 * crossprod(W[, i], R) )
      sol <- quadprog::solve.QP(Dmat = G_H, dvec = dvec, Amat = A_H, bvec = b_H, meq = 1)
      H[i, ] <- abs( sol$solution )
    }

    # Compute fitted values and SSE
    fitted <- 0.5 * ( z %*% B + W %*% H )
    sse <- sum( (x - fitted)^2 )

    # Check convergence based on SSE change
    if ( abs(sse - sse_old) < tol )  break
    sse_old <- sse
  }

  runtime <- proc.time() - runtime

  list(B = B, W = W, H = H, fitted = fitted, obj = sse, iters = iter, runtime = runtime)
}

