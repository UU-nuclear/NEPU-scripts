#' Structure Gaussian Process Model with Observations
#'
#' Constructs a Gaussian Process (GP) model to analyze and predict outcomes based on a set of
#' hyper-parameters, model configurations, observed data, and their covariance matrix. It enables
#' specifying a nugget effect for numerical stability and includes functions for updating and
#' interacting with the model.
#'
#' @param hyper_par_energy_grid Numeric vector, grid of energy values at which the amplitude hyper-parameter is defined.
#' @param hyper_pars Numeric vector, the hyper-parameters for the GP model, including amplitude (sigmas) and length scale. Should have the length length(hyper_par_energy_grid)+3, and have the structure c(sigmas, len, sigma_latent, len_latent)
#' @param model List, containing the model configuration as created using defect_model() from the defect-model.R file.
#' @param observed Numeric vector, the observed data points, normally the residual between the TALYS prediction and the experimental data.
#' @param observed_cov_mat List, holding named entries for the observed covariance matrix components:
#'        - D: a diagonal matrix with statistical uncertainties,
#'        - U: a matrix holding the systematic components,
#'        - S: a matrix that maps the systematic components in U to the data vector.
#'        The structure is created in the NEPU evaluation pipeline using the nucDataBaynet package.
#' @param nugget Numeric, the nugget effect value added to the diagonal of the covariance matrix for numerical stability. Default is 1e-04.
#'
#' @return A list of components and functions for interacting with the constructed GP model:
#'         - marginal_likelihood: Function to calculate the log-likelihood of the GP model given hyper-parameters.
#'         - get_pars_cov_mat: Function to retrieve the current parameter covariance matrix.
#'         - mapping_matrix: Matrix for mapping model parameters to observations.
#'         - conditional_mean_exp: Function to compute the conditional mean of parameters mapped to experimental energies.
#'         - conditional_covariance_exp: Function to calculate the conditional covariance matrix of parameters mapped to experimental energies.
#'         - conditional_mean_pars: Function to calculate the conditional mean of model parameters.
#'         - conditional_covariance_pars: Function to compute the conditional covariance matrix of model parameters.
#'         - mapping_matrix_pars: Matrix for mapping between model parameters and covariance matrix to match energies.
#' @examples
#' # TODO
#' @export
structureGPmodel <- function(hyper_par_energy_grid, hyper_pars, model, observed, observed_cov_mat, nugget=1e-04) {
  # observed_cov_mat should be a list holding named entries
  # D - a diagonal matrix with statistical uncertainties
  # U - a matrix holding the systematic components
  # S - a matrix that maps the systematic components in U to the data vector
  
  # define for convenience
  n_sigmas <- length(hyper_par_energy_grid)
  
  correlation_func <- function(x, length_scale) {
    d_over_l <- abs(outer(x,x,"-"))/length_scale
    
    idx <- which(d_over_l<1, arr.ind=TRUE)
    d_over_l <- d_over_l[idx]
    data <- (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)
    
    sparseMatrix(i=idx[,1], j=idx[,2], x=data)
  }
  
  make_cov_mat <- function(hyper_pars, model, nugget=1e-04) {
    sigmas <- hyper_pars[1:n_sigmas]
    len <- hyper_pars[n_sigmas + 1]
    
    model_energies <- model$parsDt[,unique(L1)]
    corr_mat_inner <-  correlation_func(model_energies, length_scale=len)
    
    sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
    sigma_scale <- outer(sigmas_tmp,sigmas_tmp)
    cov_mat <- corr_mat_inner*sigma_scale
    
    return(cov_mat + Diagonal(nrow(cov_mat),x=nugget))
  }
  
  update_cov_mat <- function(hyper_pars) {
    # if the hyper-parameters have changed from the last call
    # re-build covariance matrix and calculate its inverse
    if(!all(hyper_pars==cur_hyper_pars)) {
      cur_hyper_pars <<- hyper_pars
      cur_cov_mat <<- make_cov_mat(cur_hyper_pars, model)
      cur_cov_mat_inv <<- solve(cur_cov_mat)
      cur_cholZ <<- make_cholesky_decomp_Z()
    }
  }
  
  make_cholesky_decomp_Z <- function() {
    U_inv <- bdiag(c(list(U_obs_inv),replicate(n_channels,cur_cov_mat_inv, simplify=FALSE)))
    
    Z <- forceSymmetric(Sexp_Dinv_SexpT + U_inv)
    makeCholesky(Z)
    # Creating the cholesky decomposition seems to be the slowest part
    # It may be done faster since not the full matrix changes and the part
    # that changes is a repeated block diagonal
  }
  
  local_mult_xt_invCov_x <- function() {
    crossprod(observed, (invD_obsX - solve(observed_cov_mat$D, S_exp %*% solve(cur_cholZ, crossprod(S_exp, invD_obsX)))))
  }
  
  local_mult_invCov_x <- function() {
    invD_obsX - solve(observed_cov_mat$D, S_exp %*% solve(cur_cholZ, crossprod(S_exp, invD_obsX)))
  }
  
  local_chisquare <- function() {
    # calculate the chisquare term of the marginal log-likelihood
    as.vector(local_mult_xt_invCov_x())
  }
  
  local_logDetCov <- function() {
    tmp <- determinant(bdiag(c(list(observed_cov_mat$U), replicate(n_channels,cur_cov_mat, simplify=FALSE))))
    stopifnot(tmp$sign == 1)
    logDetU <- tmp$modulus
    
    tmp <- determinant(cur_cholZ)
    stopifnot(tmp$sign == 1)
    logDetZ <- 2 * tmp$modulus
    
    as.vector(logDetU + logDetD + logDetZ)
  }
  
  latent_cov_mat <- function(sigma, len) {
    sigma^2*correlation_func(hyper_par_energy_grid, length_scale=len)
  }
  
  log_likelihood <- function(hyper_pars) {
    
    # update the covariance matrix if hyper-parameters have changed
    update_cov_mat(hyper_pars)
    
    chi_square <- local_chisquare()
    log_det_cov <- local_logDetCov()
    L <- as.vector(-0.5*(length(observed)*log(2*pi) + log_det_cov + chi_square))
    
    # GP prior on the sigmas vector
    
    # the amplitude of the GP on sigmas covariance function
    # sigma0 <- 1 # 0.01 correspond to a 10% relative error on the talys prediction
    # the GP prior on sigmas covariance matrix:
    sigmas_prior_cov_mat <- latent_cov_mat(cur_hyper_pars[length(cur_hyper_pars)-1], cur_hyper_pars[length(cur_hyper_pars)])
    
    tmp <- determinant(sigmas_prior_cov_mat)
    stopifnot(tmp$sign == 1)
    log_det_sigmas_prior_cov_mat <- tmp$modulus
    
    sigmas <- hyper_pars[1:n_sigmas]
    chi_square <- (sigmas - sqrt(sigmas_mean)) %*% solve(sigmas_prior_cov_mat, (sigmas-sqrt(sigmas_mean)))
    
    P <- as.vector(-0.5*(length(sigmas)*log(2*pi) + log_det_sigmas_prior_cov_mat + chi_square))
    
    n_calls <<- n_calls + 1
    if((n_calls %% 25) == 0) print(paste("------- log_likelihood calls:",n_calls,":",L+P))
    
    L + P
  }
  
  get_pars_cov_mat <- function(hyper_pars) {
    # update the covariance matrix if hyper-parameters have changed
    update_cov_mat(hyper_pars)
    
    return(bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)))
  }
  
  conditional_mean_exp <- function(hyper_pars) {
    # returns the conditional mean of parameters mapped to the experimental
    # energies. 
    
    update_cov_mat(hyper_pars)
    
    # UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    # SS <- cbind(observed_cov_mat$S, mapping_matrix)
    # Sigma_yy_inv_y <- mult_invCov_x(observed, observed_cov_mat$D, SS, UU)
    # return(mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix) %*% Sigma_yy_inv_y)
    
    Sigma_yy_inv_y <- local_mult_invCov_x()
    return(mapping_matrix %*% bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)) %*% t(mapping_matrix) %*% Sigma_yy_inv_y)
  }
  
  conditional_covariance_exp <- function(hyper_pars) {
    # returns the conditional covariance matrix of parameters mapped to the
    # experimental energies. 
    update_cov_mat(hyper_pars)
    
    # UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    # SS <- cbind(observed_cov_mat$S, mapping_matrix)
    # SxKxST <- mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix)
    # return(SxKxST - mult_xt_invCov_x(SxKxST,observed_cov_mat$D, SS, UU))
    
    
    SxKxST <- mapping_matrix %*% bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)) %*% t(mapping_matrix)
    
    invDx <- solve(observed_cov_mat$D, SxKxST)
    
    return(SxKxST - crossprod(SxKxST, invDx - solve(observed_cov_mat$D, S_exp %*% solve(cur_cholZ, crossprod(S_exp,invDx)))))
  }
  
  conditional_mean_pars <- function(hyper_pars) {
    # returns the conditional mean of parameters
    update_cov_mat(hyper_pars)
    
    # UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    # SS <- cbind(observed_cov_mat$S, mapping_matrix)
    # return(mapping_matrix_pars  %*% cur_cov_mat %*% t(mapping_matrix) %*%  mult_invCov_x(observed, observed_cov_mat$D, SS, UU))
    
    Sigma_yy_inv_y <- local_mult_invCov_x()
    return(mapping_matrix_pars %*% bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)) %*% t(mapping_matrix) %*% Sigma_yy_inv_y)
  }
  
  conditional_covariance_pars <- function(hyper_pars) {
    # returns the conditional covariance matrix of parameters
    update_cov_mat(hyper_pars)
    
    # UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    # SS <- cbind(observed_cov_mat$S, mapping_matrix)
    # SxKxST <- mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix_pars)
    # return(mapping_matrix_pars %*% cur_cov_mat %*% t(mapping_matrix_pars) - mult_xt_invCov_x(SxKxST,observed_cov_mat$D, SS, UU))
    
    K_beta <- mapping_matrix_pars %*% bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)) %*% t(mapping_matrix_pars)
    JK_beta <- mapping_matrix %*% bdiag(replicate(n_channels, cur_cov_mat, simplify=FALSE)) %*% t(mapping_matrix_pars)
    invDx <- solve(observed_cov_mat$D, JK_beta)
  
    return(K_beta - crossprod(JK_beta, invDx - solve(observed_cov_mat$D, S_exp %*% solve(cur_cholZ, crossprod(S_exp, invDx)))))
    
    # A1 <- solve(cur_cholZ, crossprod(S_exp, invDx))
    # A2 <- S_exp %*% A1
    # A3 <- solve(observed_cov_mat$D, A2)
    # A4 <- invDx - A3
    # A5 <- crossprod(JK_beta, A4)
  }
  
  # create the mapping matrix
  # maps between the parameters of the model and the covariance matrix to match the energies
  # it also takes into account the scaling, to convert the value in the covariance matrix to
  # a relative variance. So multiplying this matrix (from the left) with a unit vector
  # will result in the default TALYS prediction for the exclusive channels at the energy grid
  # of the defect model (i.e. the column V1 in model$parsDt), multiplying this result with
  # model$jac will thereby result in the default TALYS prediction at the experimental points in
  # model$expDt. Sandwiching the matrix created by make_cov_mat() will result in the full
  # parameter covariance matrix K_beta
  n_channels <- length(model$parsDt[,unique(REAC)])
  model_energies <- model$parsDt[,unique(L1)]
  channel_idx = 0
  idx2 <- c()
  for(reac in model$parsDt[,unique(REAC)]) {
    target_dt <- data.table(L1 = model_energies, EnergyIndex = channel_idx + seq_along(model_energies))
    idx2 <- c(idx2, model$parsDt[REAC==reac][target_dt, on=.(L1), nomatch= 0L][,EnergyIndex])
    channel_idx <- channel_idx + length(model_energies)
  }
  idx1 <- model$parsDt[,IDX]
  value <- model$parsDt[,V1]
  
  mapping_matrix_pars <- sparseMatrix(idx1,idx2,x=value)
  
  # finally multiply together the mapping matrix and the model Jacobian
  mapping_matrix <- model$jac %*% mapping_matrix_pars
  
  S_exp <- cbind(observed_cov_mat$S, mapping_matrix)
  
  # precalculate fixed terms used in the Cholesky decomposition and the Chisquare
  invD_obsX <- solve(observed_cov_mat$D, observed)
  U_obs_inv <- solve(observed_cov_mat$U)
  Sexp_Dinv_SexpT <- crossprod(S_exp, solve(observed_cov_mat$D, S_exp))
  
  # member variables used by both log_likelihood() and grad_log_likelihood()
  cur_hyper_pars <- hyper_pars
  cur_cov_mat <- make_cov_mat(cur_hyper_pars, model)
  cur_cov_mat_inv <- solve(cur_cov_mat)
  cur_cholZ <- make_cholesky_decomp_Z()
  
  tmp <- determinant(observed_cov_mat$D)
  stopifnot(tmp$sign == 1)
  logDetD <- tmp$modulus
  
  # the mean of the GP on the sigma hyper-paramter is set to the estimate of the
  # standard deviation of the residual/default prediction, we use the more robust
  # estimate of the standard deviation obtained from 1.4826*MAD(x)
  #sigmas_mean <- 1.4826*mad(observed/(model$jac %*% model$parsDt[,V1]))
  #sigmas_mean <- 0.1
  sigmas_mean <- 0
  
  n_calls <- 0
  
  list(marginal_likelihood=log_likelihood,
       get_pars_cov_mat=get_pars_cov_mat,
       mapping_matrix = mapping_matrix,
       conditional_mean_exp = conditional_mean_exp,
       conditional_covariance_exp = conditional_covariance_exp,
       conditional_mean_pars=conditional_mean_pars,
       conditional_covariance_pars=conditional_covariance_pars,
       mapping_matrix_pars=mapping_matrix_pars)
}