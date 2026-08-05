snmfkl.eg <- function(x, k, W = NULL, H = NULL, lr_w = 0.1, lr_h = 0.1, maxiter = 1000, tol = 1e-6, clip_exp = 20) {

  runtime <- proc.time()
  n <- dim(x)[1]  ;  D <- dim(x)[2]
  # Initialize W, H randomly on simplex
  if ( is.null(W) ) {
    W <- matrix( rangen::Runif(n * k), nrow = n, ncol = k)
    W <- W / Rfast::rowsums(W)
  }
  if ( is.null(H) ) {
    H <- matrix( rangen::Runif(k * D), nrow = k, ncol = D)
    H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
  }
  prev_obj <- Inf
  Z <- W %*% H              # model

  for ( iter in 1:maxiter ) {
    R <- x / Z                # elementwise
    # ---- Gradient for KL divergence ----
    grad_w <- tcrossprod(- R, H)
    grad_h <- crossprod(W, - R)
    # ---- Exponentiated-gradient updates ----
    expo_w <-  -lr_w * grad_w
    expo_w <- pmax( pmin(expo_w, clip_exp), -clip_exp )
    W <- W * exp(expo_w)
    W <- W / Rfast::rowsums(W)
    expo_h <-  -lr_h * grad_h
    expo_h <- pmax(pmin(expo_h, clip_exp), -clip_exp)
    H <- H * exp(expo_h)
    H <- H / Rfast::rowsums(H)  ## FIX: rows sum to 1, not columns
    # ---- Compute KL divergence ----
    Z <- W %*% H
    obj <- sum( x * log(x / Z), na.rm = TRUE )
    if ( abs(prev_obj - obj ) < tol) {
      break
    }
    prev_obj <- obj
  }

  colnames(H) <- colnames(x)
  runtime <- proc.time() - runtime

  list(W = W, H = H, Z = Z, obj = obj, iters = iter, runtime = runtime)
}
