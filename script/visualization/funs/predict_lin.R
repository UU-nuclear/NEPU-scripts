predict_lin <- function(x, xp, y, cov_y, xk, mu, nugget, kernel, ...){
  
  this.mu <- if(is.function(mu)){mu(x)} else {rep(0, length(x))} 
  this.mup <- if(is.function(mu)){mu(xp)} else {rep(0, length(xp))} 
  
  d <- as.vector(y - this.mu)  
  
  S_xk <- assemble_interpolation_matrix(x, xk)
  S_xpk <- assemble_interpolation_matrix(xp, xk)
  
  K_kk <- covariance_from_kernel(x1=xk, x2=xk, kernel, ...) + amp^2 *  Diagonal(n=length(xk), x=nugget)
  k <- covariance_from_kernel(x, x, kernel, ...)
  k_xpxp <- covariance_from_kernel(xp, xp, kernel, ...)
  
  SKS_xkk <- forceSymmetric(S_xk %*% K_kk %*% t(S_xk)) 
  D <- diag(diag(k) - diag(SKS_xkk)) + cov_y + amp^2 * Diagonal(n=length(y), x=nugget) 
  covi <- woodburry(D, S_xk, K_kk)
  
  S_xpx <- S_xpk %*% K_kk %*% t(S_xk)
  
  mu_xp <- this.mup + S_xpx %*% covi %*% d
  var_fxp <- k_xpxp - S_xpx %*% covi %*% t(S_xpx)
  
  
  list("mu"=mu_xp,
       "var"=var_fxp)
}