#  nucdataBaynet - Nuclear Data Evaluation Using a Bayesian Network 
#  Copyright (C) 2019  Georg Schnabel
#  
#  nucdataBaynet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  nucdataBaynet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>


#' Levenberg-Marquardt algorithm with prior
#' extended with parallel function evaluation
#'
#' @param fn model function
#' @param jac Jacobian of model function
#' @param pinit initial model parameters
#' @param p0 prior mean of model parameters
#' @param P0 prior covariance matrix of model parameters
#' @param yexp experimental data points
#' @param D only diagonal contribution to experimental covariance matrix
#' @param S mapping from systematic components to observations
#' @param X prior covariance matrix for systematic components
#' @param control control arguments, see \code{details}
#'
#' @details
#' The control argument is a list that can contain
#' \itemize{
#'   \item maximum number of iterations \code{maxit}
#'   \item absolute tolerance required for convergence \code{abstol}
#'   \item relative tolerance required for convergence \code{reltol}
#' }
#'
#' @return
#' list with elements \code{par}, \code{fn}, \code{jac}, \code{counts}, \code{value}
#' @export
#'
LMalgo_parallel <- function(fn, jac, pinit, p0, P0, yexp, D, S, X,
                   lower = -Inf, upper = Inf, logger = NULL,
                   control = list(maxit = 10, abstol = 1e-3,
                                  reltol = 1e-3, steptol = 1e-4, mu = NULL,
                                  nproc = 1, strategy = "likelihood"),
                   J = NULL
                   ) {

  parDefaults <- list(
    control = list(maxit = 10, abstol = 1e-3, reltol = 1e-3)
  )
  control <- modifyList(parDefaults$control, control)

  isDebugOn <- TRUE
  debuginfo <- function(str) {
    if (isDebugOn) cat(str, "\n")
  }

  # initialization
  nproc <- max(3,control$nproc)
  tau <- 1e-5  # factor to determine damping factor
  DD <- Diagonal(x = rep(1, length(pinit)))

  # speeds up inversions of the experimental covariance matrix
  cholZ <- makeCholZ(D, S, X)
  breakCounter <- 0

  invP0 <- solve(P0)
  pref <- as.vector(pinit)

  if(is.null(J)) {
    debuginfo("calculating initial jacobian...")
    J <- jac(as.vector(pref))
  }
  fref <- fn(as.vector(pref))

  dpriorRef <- pref - p0
  LpriorRef <- as.vector(crossprod(dpriorRef, invP0 %*% dpriorRef))
  Lref <- as.vector(mult_xt_invCov_x(yexp - fref, D, S, X, cholZ = cholZ)) + LpriorRef

  logBuffer <- list()
  i <- 0  # counter
  accepted <- TRUE # checking if the last iteration lead to an accepted step

  tJinvBJ <- forceSymmetric(mult_xt_invCov_x(J, D, S, X, cholZ = cholZ))
  invP1 <- invP0 + tJinvBJ
  mu <- if (is.null(control$mu)) tau * max(diag(invP1)) else control$mu

  while (i < control$maxit && breakCounter < 2) {
    i <- i + 1
    debuginfo(paste0("iteration: ",i))

    # propose new parameter set and calculate function value
    tJinvBJ <- forceSymmetric(mult_xt_invCov_x(J, D, S, X, cholZ = cholZ))
    invP1 <- invP0 + tJinvBJ

    n1 <- floor(0.5*nproc)
    n2 <- ceiling(0.5*nproc) + 2
    mus <- c(seq(n1,1),1/seq(3,n2))*mu

    if(!accepted) {
      # try only larger values of mu if the former iteration did not find an improvement 
      mus <- seq(nproc,1)*mu
    }

    #prepare the next set of proposal parameters
    pprops <- matrix(0,ncol=length(mus),nrow=length(pref))
    for(col in 1:length(mus)) {
          m <- t(J) %*% mult_invCov_x(yexp - fref, D, S, X, cholZ = cholZ) + invP0 %*% (p0 - pref)
          pprop <- as.vector(pref + solve(invP1 + mus[col] * DD, m))

           # project on box constraints
          pprop_unconstr <- pprop
          pprop[pprop < lower] <- lower[pprop < lower]
          pprop[pprop > upper] <- upper[pprop > upper]

          pprops[,col] <- pprop
    }
   
    debuginfo("calculating function values...")
    # calculate the function values at the proposal parameter sets
    fprops <- fn(pprops)
    debuginfo("...done!")
    
    Lprop_min <- Inf
    Lprop_approx_min <- Inf
    fprop_approx_min <- Inf
    col_min <- 0

    former_step_gain <- 0
    for(col in length(mus):1) {
      fprop <- fprops[,col]
      pprop <- pprops[,col]
      fprop_approx <- fref + J %*% (pprop - pref)
      fprop_approx_unconstr <- fref + J %*% (pprop_unconstr - pref)

      # calculate true objective function and approximation thereof
      dpriorProp <- pprop - p0
      dpriorProp_unconstr <- pprop_unconstr - p0

      LpriorProp <- as.vector(crossprod(dpriorProp, invP0 %*% dpriorProp))
      LpriorProp_unconstr <- as.vector(crossprod(dpriorProp_unconstr, invP0 %*% dpriorProp_unconstr))

      Lprop <- as.vector(mult_xt_invCov_x(yexp - fprop, D, S, X, cholZ = cholZ)) + LpriorProp
      Lprop_approx <- as.vector(mult_xt_invCov_x(yexp - fprop_approx, D, S, X, cholZ = cholZ)) + LpriorProp
      Lprop_approx_unconstr <- as.vector(mult_xt_invCov_x(yexp - fprop_approx_unconstr, D, S, X, cholZ = cholZ)) + LpriorProp_unconstr

      #print(paste(mus[col],Lprop,Lprop_approx,Lref))
      if(control$strategy=="likelihood") {
        # this strategy picks the step with the largest local decrease in L.
        # the strategy is relatively fast, but has the risk of slowing down close to the minimum
        # see [Marquadt, Jour. Soc. for Industrial and Applied Mathematics, vol. 11, pp. 431441, June 1963]
        # if multiple minima occurs in every step picking the deepest one may lead to unwanted behaviour
        if(Lprop < Lprop_min) {
          Lprop_min <- Lprop
          Lprop_approx_min <- Lprop_approx
          fprop_approx_min <-fprop_approx
          col_min <- col
          mu <- mus[col]
        }
      }
      if(control$strategy=="gain") {
        # this strategy picks the step with the "best" gain. Starting from the smallest values
        # of mu it goes in order of increasing mu-value, and stops once the mu-value goes
        # above a threshold value of 0.75. If all values are above 0.75 the one with the smallest
        # mu is taken. If all values are below 0.75 the largest value of mu is taken
        # this will lead to a maximum step length while still in the region where the 
        # linear approximation is acceptable. If the gain is close to 1, we are well inside this region.
        # As gain -> 0 we are moving outside of this region.
        
        step_gain <- ((Lref - Lprop)+1e-10) / (abs(Lref - Lprop_approx)+1e-10)
        #print(paste(col,mu,step_gain))
        if(step_gain>0.75) {
          if((step_gain-0.75) < (0.75 - former_step_gain)) {
            # if current step gain is closer to 0.75 than the former one
            Lprop_min <- Lprop
            Lprop_approx_min <- Lprop_approx
            fprop_approx_min <-fprop_approx
            col_min <- col
            mu <- mus[col]
            break
          } else {
            # otherwise, the former_step_gain is closer to 0.75
            break
          }
        }

        Lprop_min <- Lprop
        Lprop_approx_min <- Lprop_approx
        fprop_approx_min <-fprop_approx
        col_min <- col
        mu <- mus[col]
      }
      
    }

    # +1e-10 to avoid problems if Lprop_approx == Lref
    # abs in denominator important because of
    # projecting back into feasible region
    # and then Lprop_approx < Lref not guaranteed anymore

    gain <- ((Lref - Lprop_min)+1e-10) / (abs(Lref - Lprop_approx_min)+1e-10)
    #print(paste("gain = ",gain))

    stepLength <- norm(as.matrix(pprops[,col_min]-pref),"f")
    longestStep <- max(abs(pprops[,col_min]-pref))

    # check break conditions
    #if (abs(Lprop_min - Lref) / abs(Lref) < control$reltol ||
    #    abs(Lprop_min - Lref) < control$abstol ||
    #    longestStep < control$steptol) {
    if(longestStep < control$steptol) {
      #print("   break condition TRUE")
      breakCounter <- breakCounter + 1
    } else {
      #print("   break condition FALSE")
      breakCounter <- 0
    }


    # print status information
    logBuffer <- list(iteration = i, mu = mu, gain = gain, pref = pref, fref = fref, Lref = Lref,
                      pprop = pprops[,col_min], fprop = fprops[,col_min], Lprop = Lprop_min, 
                      fprop_approx = fprop_approx_min, Lprop_approx = Lprop_approx_min,
                      Jref = J, p0 = p0, P0 = P0, yexp = yexp, D = D, S = S, X = X, 
                      stepLength = stepLength, longestStep = longestStep)
    
    if (is.function(logger)) logger(logBuffer)
    
    # accept if proposed parameter set better than old one
    accepted <- FALSE
    if (Lprop_min < Lref) {
      pref <- pprops[,col_min]
      fref <- fprops[,col_min]

      debuginfo("calculating Jacobian...")
      J <- jac(as.vector(pprops[,col_min]))
      debuginfo("...done!")
      
      Lref <- Lprop_min
      accepted <- TRUE
    } else {
      # set the next mu-value to 2*(current mu-value)
      mu <- 2*mu
    }
  }
  
  tJinvBJ <- forceSymmetric(mult_xt_invCov_x(J, D, S, X, cholZ = cholZ))
  invP1 <- invP0 + tJinvBJ
  
  parCov <- solve(invP1)
  stdAsyFitErr <- sqrt(diag(J %*% parCov %*% t(J)))
  
  list(par = pref,
       fn = fref,
       jac = J,
       counts = i,
       value = Lref,
       parCovLM = parCov,
       stdAsyFitErrLM = stdAsyFitErr
  )
}









