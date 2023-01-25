# transformations used to restrict TALYS parameters 
# to lie in a certain interval. The logistic function
# is used as transformation.
# Let xt = f(x) with the transformed parameter xt and the original parameter x 
# Parameters:
#  x0    ... the turning point where f(x0) = x0, i.e. transformed parameter 
#            equals untransformed parameter 
#  x_min ... is the lower limit for the original parameters 
#  x_max ... is the upper limit for the original parameters 
# the rate of change is adapted so that df/dx = 1 at x = x0
#
# parameterTransform expects vectors x0, x_min, x_max, each of the same
# length as the number of parameters holding the limits of the corresponding
# parameter
# 
# erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
# 
# erfinv <- function (x) qnorm((1 + x)/2)/sqrt(2)
#

transFunc <- function(x, x0, delta) {
  x0 + delta*(2*pnorm((x-x0)/delta) - 1)
}

inv_transFunc <- function(z, x0, delta) {
  x0 + delta*qnorm(0.5*(1+(z-x0)/delta))
}

jac_transFunc <- function(x, x0, delta) {
  x <- 2*delta*dnorm(x,x0,delta)
  diag(x = x, nrow = length(x))
}

parameterTransform <- function(x0, delta) {

  fun <- function(x) {
    if(is.vector(x)) {
      transFunc(x, x0, delta)
    } else if(is.matrix(x)) {
      apply(x,2,transFunc,x0=x0,delta)
    }
  }

  invfun <- function(z) {
    if(is.vector(z)) {
      inv_transFunc(z, x0, delta)
    } else if(is.matrix(z)) {
      apply(z,2,inv_transFunc,x0=x0,delta)
    }
  }

  jac <- function(x) {
    jac_transFunc(x, x0, delta)
  }

  list(fun = fun, jac = jac, invfun = invfun)
}
