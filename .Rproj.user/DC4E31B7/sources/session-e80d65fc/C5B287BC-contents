# Define helper OUTSIDE the main function to avoid closure capture
.solve_w_row_qp <- function(i, G_W, g_W) {
  nnsolve::fnnls(G_W, g_W[, i], sum_to_constant = TRUE, constant = 1)
}

snmf.qp <- function(x, k, W = NULL, H = NULL, k_means = TRUE, lr_h = 0.1,
                    maxiter = 1000, tol = 1e-6, ridge = 1e-8, history = FALSE, ncores = 1) {

  runtime <- proc.time()
  n <- dim(x)[1]  ;  D <- dim(x)[2]
  if ( n <= D )  veo <- TRUE

  # Initialize W, H on simplex
  if ( k_means ) {
    H <- kmeans(x, k)$centers
    if ( veo ) {
      H <- matrix(nrow = k, ncol = D)
      cs <- as.numeric( sparcl::KMeansSparseCluster(x, k, wbounds = 2)[[1]]$Cs )
      for ( i in 1:k ) H[i, ] <- Rfast::colmeans(x[cs == i, , drop = FALSE])
    }
  } else {
    if ( is.null(H) ) {
      H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D )
      H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
    }
  }

  W <- matrix(nrow = n, ncol = k)
  error <- numeric(maxiter)
  ridgek <- diag(ridge, k)

  suppressWarnings({

    # CREATE CLUSTER ONCE BEFORE LOOP (like snmf.sqp)
    if (ncores > 1) {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      parallel::clusterEvalQ(cl, library(nnsolve))
      parallel::clusterExport( cl, varlist = c(".solve_w_row_qp", "G_W", "g_W"), envir = environment() )
    }

    if ( !veo ) {  ## veo is FALSE, n > p

      sx2 <- sum(x^2)
      kD <- k * D
      XX <- matrix(0, kD, kD)
      ind <- matrix(1:kD, ncol = k, byrow = TRUE)
      A <- XX
      for (i in 1:k)  A[i, ind[, i]] <- 1
      A <- t( rbind(A, diag(kD)) )
      A <- A[, -c((k + 1):kD)]
      bvec <- c( rep(1, k), rep(0, kD) )

      for ( it in 1:maxiter ) {
        # ----------------- W update -----------------
        G_W <- 2 * tcrossprod(H) + ridgek
        g_W <- 2 * tcrossprod(H, x)

        if ( ncores > 1 ) {
          # Export iteration-specific variables
          parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment())
          W <- t( parallel::parSapply(cl, 1:n, .solve_w_row_qp, G_W = G_W, g_W = g_W) )
        } else {
          for ( i in 1:n )  W[i, ] <- nnsolve::fnnls(G_W, g_W[, i], sum_to_constant = TRUE, constant = 1)
        }

        # ----------------- H update -----------------
        dvec <- as.vector( crossprod(W, x) )
        xx <- crossprod(W)
        for ( i in 1:D )  XX[ind[i, ], ind[i, ]] <- xx
        f <- try( quadprog::solve.QP(Dmat = XX, dvec = dvec, Amat = A, bvec = bvec, meq = k), silent = TRUE )
        if ( identical(class(f), "try-error") ) {
          f <- quadprog::solve.QP(Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = k)
        }
        H <- matrix(abs(f$solution), ncol = D)
        err <- sx2 + 2 * f$value
        error[it] <- err

        if ( it > 1 && abs(error[it - 1] - err) < tol ) {
          break
        }
      }
      Z <- W %*% H

    } else {  ## veo is TRUE, n < p

      for ( it in 1:maxiter ) {
        # ----------------- W update -----------------
        G_W <- 2 * tcrossprod(H) + ridgek
        g_W <- 2 * tcrossprod(H, x)

        if ( ncores > 1 ) {
          # Export iteration-specific variables
          parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment())
          W <- t( parallel::parSapply(cl, 1:n, .solve_w_row_qp, G_W = G_W, g_W = g_W) )
        } else {
          for ( i in 1:n )  W[i, ] <- nnsolve::fnnls(G_W, g_W[, i], tol = tol, sum_to_constant = TRUE, constant = 1)
        }

        # ----------------- H update -----------------
        E <- W %*% H - x
        grad_h <- crossprod(W, E)
        H <- H * exp(-lr_h * grad_h)
        H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
        Z <- W %*% H
        err <- sum( (x - Z)^2 )
        error[it] <- err

        if ( it > 1 && abs(error[it - 1] - err) < tol ) {
          break
        }
      }

    }  ##  end if (!veo)

  })  # end suppressWarnings

  runtime <- proc.time() - runtime
  error <- error[1:it]
  obj <- error[it]
  if ( !history )  error <- NULL
  colnames(H) <- colnames(x)

  list(W = W, H = H, Z = Z, obj = obj, error = error, iters = it, runtime = runtime)
}



# # Define helper OUTSIDE the main function to avoid closure capture
# .solve_w_row_qp <- function(i, G_W, g_W, A_W, b_W) {
#   sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
#   abs(sol$solution)
# }
#
# snmf.qp_old <- function(x, k, W = NULL, H = NULL, k_means = TRUE, veo = FALSE, lr_h = 0.1,
#                      maxiter = 1000, tol = 1e-6, ridge = 1e-8, history = FALSE, ncores = 1) {
#
#   runtime <- proc.time()
#   n <- dim(x)[1]  ;  D <- dim(x)[2]
#
#   # Initialize W, H on simplex
#   if ( k_means ) {
#     H <- kmeans(x, k)$centers
#     if ( veo ) {
#       H <- matrix(nrow = k, ncol = D)
#       cs <- as.numeric( sparcl::KMeansSparseCluster(x, k, wbounds = 2)[[1]]$Cs )
#       for ( i in 1:k ) H[i, ] <- Rfast::colmeans(x[cs == i, , drop = FALSE])
#     }
#   } else {
#     if ( is.null(H) ) {
#       H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D )
#       H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
#     }
#   }
#
#   W <- matrix(nrow = n, ncol = k)
#   error <- numeric(maxiter)
#   A_W <- t( rbind(rep(1, k), diag(k)) )
#   b_W <- c(1, rep(0, k))
#   ridgek <- diag(ridge, k)
#
#   suppressWarnings({
#
#   # CREATE CLUSTER ONCE BEFORE LOOP (like snmf.sqp)
#   if (ncores > 1) {
#     cl <- parallel::makeCluster(ncores)
#     on.exit(parallel::stopCluster(cl), add = TRUE)
#     parallel::clusterEvalQ(cl, library(quadprog))
#     parallel::clusterExport( cl, varlist = c(".solve_w_row_qp", "A_W", "b_W"), envir = environment() )
#   }
#
#   if ( !veo ) {  ## veo is FALSE, n > p
#
#     sx2 <- sum(x^2)
#     kD <- k * D
#     XX <- matrix(0, kD, kD)
#     ind <- matrix(1:kD, ncol = k, byrow = TRUE)
#     A <- XX
#     for (i in 1:k)  A[i, ind[, i]] <- 1
#     A <- t( rbind(A, diag(kD)) )
#     A <- A[, -c((k + 1):kD)]
#     bvec <- c( rep(1, k), rep(0, kD) )
#
#     for ( it in 1:maxiter ) {
#       # ----------------- W update -----------------
#       G_W <- 2 * tcrossprod(H) + ridgek
#       g_W <- 2 * tcrossprod(H, x)
#
#       if ( ncores > 1 ) {
#         # Export iteration-specific variables
#         parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment())
#         W <- t( parallel::parSapply(cl, 1:n, .solve_w_row_qp, G_W = G_W, g_W = g_W, A_W = A_W, b_W = b_W) )
#       } else {
#         for ( i in 1:n ) {
#           sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
#           W[i, ] <- abs(sol$solution)
#         }
#       }
#
#       # ----------------- H update -----------------
#       dvec <- as.vector( crossprod(W, x) )
#       xx <- crossprod(W)
#       for ( i in 1:D )  XX[ind[i, ], ind[i, ]] <- xx
#       f <- try( quadprog::solve.QP(Dmat = XX, dvec = dvec, Amat = A, bvec = bvec, meq = k), silent = TRUE )
#       if ( identical(class(f), "try-error") ) {
#         f <- quadprog::solve.QP(Dmat = Matrix::nearPD(XX)$mat, dvec = dvec, Amat = A, bvec = bvec, meq = k)
#       }
#       H <- matrix(abs(f$solution), ncol = D)
#       err <- sx2 + 2 * f$value
#       error[it] <- err
#
#       if ( it > 1 && abs(error[it - 1] - err) < tol ) {
#         break
#       }
#     }
# 	Z <- W %*% H
#
#   } else {  ## veo is TRUE, n < p
#
#     sx2 <- sum(x^2)
#     kD <- k * D
#     XX <- matrix(0, kD, kD)
#     ind <- matrix(1:kD, ncol = k, byrow = TRUE)
#     A <- XX
#     for (i in 1:k)  A[i, ind[, i]] <- 1
#     A <- t( rbind(A, diag(kD)) )
#     A <- A[, -c((k + 1):kD)]
#     bvec <- c( rep(1, k), rep(0, kD) )
#
#     for ( it in 1:maxiter ) {
#       # ----------------- W update -----------------
#       G_W <- 2 * tcrossprod(H) + ridgek
#       g_W <- 2 * tcrossprod(H, x)
#
#       if ( ncores > 1 ) {
#         # Export iteration-specific variables
#         parallel::clusterExport(cl, varlist = c("G_W", "g_W"), envir = environment())
#         W <- t( parallel::parSapply(cl, 1:n, .solve_w_row_qp, G_W = G_W, g_W = g_W, A_W = A_W, b_W = b_W) )
#       } else {
#         for ( i in 1:n ) {
#           sol <- quadprog::solve.QP(Dmat = G_W, dvec = g_W[, i], Amat = A_W, bvec = b_W, meq = 1)
#           W[i, ] <- abs(sol$solution)
#         }
#       }
#
#       # ----------------- H update -----------------
#       E <- W %*% H - x
#       grad_h <- crossprod(W, E)
#       H <- H * exp(-lr_h * grad_h)
#       H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
#       Z <- W %*% H
#       err <- sum( (x - Z)^2 )
#       error[it] <- err
#
#       if ( it > 1 && abs(error[it - 1] - err) < tol ) {
#         break
#       }
#     }
#
#   }  ##  end if (!veo)
#
#   })  # end suppressWarnings
#
#   runtime <- proc.time() - runtime
#   error <- error[1:it]
#   obj <- error[it]
#   if ( !history )  error <- NULL
#   colnames(H) <- colnames(x)
#
#   list(W = W, H = H, Z = Z, obj = obj, error = error, iters = it, runtime = runtime)
# }



