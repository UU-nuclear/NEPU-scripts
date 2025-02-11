
GP_model <- function(hyper_pars, model, observed, observed_cov_mat, nugget=1e-04) {

	cov_mat_derivatives <- function(model, sigma, len) {

		dcov_func_dsigma <- function(x1, x2) {
			# note that this is independent of the value of sigma!!!
			d_over_l <- abs(x1-x2)/len
			ifelse(d_over_l < 1, (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi), 0)
		}

		dcov_func_dl <- function(x1, x2) {
			d_over_l <- abs(x1-x2)/len
			ifelse(d_over_l<1, (4/3)*sigma*(pi*(1-d_over_l)*cos(pi*d_over_l) + sin(pi*d_over_l))*sin(pi*d_over_l)*(d_over_l/len), 0)
		}


		dCovMat_dsigma <- list()
		dCovMat_dlen <- list()
		for(reac in model$parsDt[,unique(REAC)]) {
			model_energies <- model$parsDt[REAC==reac,L1]
			model_default <- model$parsDt[REAC==reac,V1]
			cov_mat_dsigma <- Matrix(outer(model_energies,model_energies,dcov_func_dsigma))
			cov_mat_dlen <- Matrix(outer(model_energies,model_energies,dcov_func_dl))

			scale <- outer(model_default,model_default)

			dCovMat_dsigma <- append(dCovMat_dsigma,list(Matrix(cov_mat_dsigma*scale)))
			dCovMat_dlen <- append(dCovMat_dlen,list(Matrix(cov_mat_dlen*scale)))
		}

		full_dCovMat_dsigma <- bdiag(dCovMat_dsigma)
		full_dCovMat_dlen <- bdiag(dCovMat_dlen)

		exp_dCovMat_dsigma <- model$jac %*% full_dCovMat_dsigma %*% t(model$jac)
		exp_dCovMat_dlen <- model$jac %*% full_dCovMat_dlen %*% t(model$jac)

		list(sigma=exp_dCovMat_dsigma, len=exp_dCovMat_dlen)
	}

	Trace <- function(A,B) {
		sum(A*B)
	}

	make_cov_mat <- function(hyper_pars, model, nugget=1e-04) {

		sigma <- hyper_pars[1]
		len <- hyper_pars[2]
		cov_func <- function(x1, x2) {
			d_over_l <- abs(x1-x2)/len
			ifelse(d_over_l < 1, sigma * ( (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)), 0)
		}

		priors <- list()
		for(reac in model$parsDt[,unique(REAC)]) {
			model_energies <- model$parsDt[REAC==reac,L1]
			model_default <- model$parsDt[REAC==reac,V1]
			cov_mat <- Matrix(outer(model_energies,model_energies,cov_func))

			scale <- outer(model_default,model_default)

			cov_mat <- Matrix(cov_mat*scale)

			priors <- append(priors,list(cov_mat))
		}

		full_prior_cov_mat <- bdiag(priors)
		full_prior_cov_mat + Diagonal(n=nrow(full_prior_cov_mat), x=nugget)
	}

	update_cov_mat <- function(hyper_pars) {

		# print("update_cov_mat!")
		# if the hyper-parameters have changed from the last call
		# re-build covariance matrix
		if(!all(hyper_pars==cur_hyper_pars)) {
			# print("hyper-pars have changed!")
			cur_hyper_pars <<- hyper_pars
			cur_cov_mat <<- make_cov_mat(cur_hyper_pars, model)
			cur_cov_mat_exp <<-  model$jac %*% cur_cov_mat %*% t(model$jac) + observed_cov_mat
		}

	}

	log_likelihood <- function(hyper_pars) {

		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		chi_square <- observed %*% solve(cur_cov_mat_exp, observed)

		tmp <- determinant(cur_cov_mat_exp)
		stopifnot(tmp$sign == 1)
		log_det_cov <- tmp$modulus

		as.vector(-0.5*(length(observed)*log(2*pi) + log_det_cov + chi_square))
	}

	grad_log_likelihood <- function(hyper_pars) {

		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		inv_sigma_obs <- solve(cur_cov_mat_exp)
		inv_sigma_obs_xt <- crossprod(inv_sigma_obs, observed)
		
		derivatives <- cov_mat_derivatives(model, hyper_pars[1], hyper_pars[2])

		dL_dsigma <- as.vector(-0.5 * (Trace(inv_sigma_obs,derivatives$sigma) - t(inv_sigma_obs_xt) %*% derivatives$sigma %*% inv_sigma_obs_xt))
		dL_dl <- as.vector(-0.5 * (Trace(inv_sigma_obs,derivatives$len) - t(inv_sigma_obs_xt) %*% derivatives$len %*% inv_sigma_obs_xt))

		#rbind(dL_dsigma,dL_dl)
		c(dL_dsigma,dL_dl)
	}

	get_pars_cov_mat <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		return(cur_cov_mat)
	}

	get_exp_cov_mat <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		return(cur_cov_mat_exp)
	}

	# member variables used by both log_likelihood() and grad_log_likelihood()
	cur_hyper_pars <- hyper_pars
	cur_cov_mat <- make_cov_mat(cur_hyper_pars, model)
	cur_cov_mat_exp <-  model$jac %*% cur_cov_mat %*% t(model$jac) + observed_cov_mat

	list(marginal_likelihood=log_likelihood, grad_ML=grad_log_likelihood, get_pars_cov_mat=get_pars_cov_mat, get_exp_cov_mat=get_exp_cov_mat)
}