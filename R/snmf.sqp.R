# Define helper OUTSIDE the main function to avoid closure capture
.solve_w_row <- function(i, G_W, g_W, A_W, l_W, u_W, osqp_settings) {
  sol <- osqp::solve_osqp( P = G_W, q = -g_W[, i], A = A_W,
               l = l_W, u = u_W, pars = osqp_settings )
  pmax(sol$x, 0)
}

snmf.sqp <- function(x, k, W = NULL, H = NULL, maxiter = 1000,
                     tol = 1e-4, ridge = 1e-8, history = FALSE, ncores = 1) {
  runtime <- proc.time()
  n <- dim(x)[1]  ;  D <- dim(x)[2]

  # Initialize W, H on simplex
  H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D)
  H <- H / Rfast::rowsums(H)
  W <- matrix(nrow = n, ncol = k)

  error <- numeric(maxiter)
  ridgek <- diag(ridge, k)
  sx2 <- sum( Matrix::colSums(x^2) )
  kD <- k * D
  ind <- matrix( 1:kD, ncol = k, byrow = TRUE )

  # Constraint matrix for W: sum(w) = 1, w >= 0
  A_W <- rbind( rep(1, k), diag(k) )
  l_W <- c(1, rep(0, k) )
  u_W <- c(1, rep(1, k) )

  # Constraint matrix for H
  i_sum <- rep(1:k, times = D)
  j_sum <- as.vector(ind[ind > 0])
  sum_constraints <- Matrix::sparseMatrix( i = i_sum, j = j_sum, x = rep(1, length(i_sum)), dims = c(k, kD) )
  A_H <- rbind( sum_constraints, Matrix::Diagonal(kD) )
  l_H <- c( rep(1, k), rep(0, kD) )
  u_H <- c( rep(1, k), rep(1, kD) )

  # OSQP settings
  osqp_settings <- osqp::osqpSettings(verbose = FALSE, eps_abs = 1e-8, eps_rel = 1e-8)

  # CREATE CLUSTER ONCE BEFORE LOOP
  if (ncores > 1) {
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, library(osqp))
    parallel::clusterExport(cl,
    varlist = c("A_W", "l_W", "u_W", "osqp_settings", ".solve_w_row"), envir = environment() )
  }

  for ( it in 1:maxiter ) {
    # ----------------- W update -----------------
    G_W <- as.matrix( 2 * tcrossprod(H) + ridgek )
    g_W <- as.matrix( 2 * tcrossprod(H, x) )

    if ( ncores > 1 ) {
      parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment() )
      # Use external helper function - no closure capture!
      suppressWarnings({
        W <- t( parallel::parSapply(cl, 1:n, .solve_w_row, G_W = G_W, g_W = g_W,
                A_W = A_W, l_W = l_W, u_W = u_W, osqp_settings = osqp_settings) )
      })
    } else {
      for ( i in 1:n ) {
        sol <- osqp::solve_osqp( P = G_W, q = -g_W[, i], A = A_W, l = l_W, u = u_W, pars = osqp_settings )
        W[i, ] <- pmax(sol$x, 0)
      }
    }

    # ----------------- H update -----------------
    dvec <- as.vector( crossprod(W, x) )
    xx <- crossprod(W) + ridgek
    XX <- Matrix::bdiag( replicate(D, xx, simplify = FALSE) )
    f <- osqp::solve_osqp( P = XX, q = -dvec, A = A_H, l = l_H, u = u_H, pars = osqp_settings )
    H <- matrix( pmax(f$x, 0), ncol = D)
    # Compute error
    err <- sx2 - 2 * sum(dvec * f$x) + sum(H * (xx %*% H))
    error[it] <- err

    if ( it > 1 && abs(error[it - 1] - err) < tol ) {
      break
    }
  }

  error <- error[1:it]
  obj <- error[it]
  if ( !history )  error <- NULL
  colnames(H) <- colnames(x)
  runtime <- proc.time() - runtime

  list( W = W, H = H, obj = obj, error = error, iters = it, runtime = runtime )
}

