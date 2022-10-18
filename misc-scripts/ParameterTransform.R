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

transFunc <- function(x, x0, x_min, x_max) {
  much_lower <- -2.*(x-x0)/(x0-x_min) > 1e2
  much_higher <- -2.*(x-x0)/(x_max-x0) < -1e2
  x[much_lower] <- x_min[much_lower]
  x[much_higher] <- x_max[much_higher]
  

  lower <- x <= x0 & !much_lower
  upper <- x >  x0 & !much_higher
  x[lower] <- x0[lower] + (x0[lower]-x_min[lower])*(1.-exp(-2.*(x[lower]-x0[lower])/(x0[lower]-x_min[lower])))/(1.+exp(-2.*(x[lower]-x0[lower])/(x0[lower]-x_min[lower])))
  x[upper] <- x0[upper] + (x_max[upper]-x0[upper])*(1.-exp(-2.*(x[upper]-x0[upper])/(x_max[upper]-x0[upper])))/(1.+exp(-2.*(x[upper]-x0[upper])/(x_max[upper]-x0[upper])))
  x
}

inv_transFunc <- function(z, x0, x_min, x_max) {
  z[z <= x_min] <- x_min[z < x_min] + 1e-09 # can happen due to numerical inaccuracy
  z[z >= x_max] <- x_max[z > x_max] - 1e-09 # can happen due to numerical inaccuracy

  lower <- z <= x0
  upper <- z > x0
  z[lower] <- x0[lower] + (x0[lower]-x_min[lower])/2 * log((x_min[lower] - z[lower]) / (z[lower] - 2*x0[lower] + x_min[lower]))
  z[upper] <- x0[upper] + (x_max[upper]-x0[upper])/2 * log((2*x0[upper] - x_max[upper] - z[upper]) / (z[upper] - x_max[upper]))
  z
}

jac_transFunc <- function(x, x0, x_min, x_max) {
  much_lower <- -2.*(x-x0)/(x0-x_min) > 1e2
  much_higher <- -2.*(x-x0)/(x_max-x0) < -1e2
  x[much_lower] <- double.xmin
  x[much_higher] <- double.xmin

  lower <- x <= x0
  upper <- x >  x0
  c = exp(-2*(x[lower]-x0[lower])/(x0[lower]-x_min[lower]))
  x[lower] <- 2*c^2/(1.+c)^2
  c = exp(-2*(x[upper]-x0[upper])/(x_max[upper]-x0[upper]))
  x[upper] <- 2*c^2/(1.+c)^2
  diag(x = x, nrow = length(x))
}

parameterTransform <- function(x0, x_min, x_max) {

  fun <- function(x) {
    transFunc(x, x0, x_min, x_max)
  }

  invfun <- function(z) {
    inv_transFunc(z, x0, x_min, x_max)
  }

  jac <- function(x) {
    jac_transFunc(x, x0, x_min, x_max)
  }

  list(fun = fun, jac = jac, invfun = invfun)
}

