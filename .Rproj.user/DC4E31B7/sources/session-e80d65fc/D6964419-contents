snmf.eg <- function(x, k, W = NULL, H = NULL, lr_w = 0.1, lr_h = 0.1, maxiter = 10000, tol = 1e-6) {
  # x: n x D compositional data (rows sum to 1)
  # k: number of archetypes
  # lr_w, lr_h: EG learning rates for W and H
  # returns W (n x k), H (k x D), Z = W %*% H, objective trace info
  runtime <- proc.time()
  n <- dim(x)[1]  ;   D <- dim(x)[2]
  # initialize W and H randomly and normalize
  if ( is.null(W) ) {
    W <- matrix( rangen::Runif(n * k), nrow = n, ncol = k)
    W <- W / Rfast::rowsums(W)
  }
  if ( is.null(H) ) {
    H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D)
    H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
  }
  prev_obj <- Inf
  obj <- Inf
  Z <- W %*% H           # n x D

  for ( iter in 1:maxiter ) {
    E <- Z - x             # n x D
    grad_w <- tcrossprod(E, H)
    W <- W * exp(- lr_w * grad_w)
    W <- W / Rfast::rowsums(W) # +eps?
    Z <- W %*% H
    E <- Z - x
    grad_h <- crossprod(W, E)
    H <- H * exp(- lr_h * grad_h)
    H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
    Z <- W %*% H
    obj <- sum( (x - Z)^2 )
    # convergence check
    if ( abs(prev_obj - obj) < tol ) {
      break
    }
    prev_obj <- obj
  }

  colnames(H) <- colnames(x)
  runtime <- proc.time() - runtime

  list(W = W, H = H, Z = Z, obj = obj, iters = iter, runtime = runtime)
}
