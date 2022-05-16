# transformations used to restrict TALYS parameters 
# to lie in a certain interval. The logistic function
# is used as transformation.
# Let xt = f(x) with the transformed parameter xt and the original parameter x 
# Parameters:
#  x0    ... the turning point where f(x0) = x0, i.e. transformed parameter 
#            equals untransformed parameter 
#  delta ... (x0-delta,x0+delta) is the range for the original parameters 
#  k     ... the speed at which the transition from x0-delta to x0+delta
#            around x0 happens
#  If df/dx(x0) = 1 is desired, then the condition delta*k = 2 must hold

trafoFun <- function(x, x0, delta, k) {

  x0 + delta * (1 - exp(-k*(x-x0))) / (1 + exp(-k*(x-x0)))
}

invTrafoFun <- function(z, x0, delta, k) {

  x0 + 1/k * log((x0 - delta - z) / (z - x0 - delta))
}

trafoJac <- function(x, x0, delta, k) {

  s <- (x-x0)
  c <- exp(-k*s)
  diag(x = delta * (k*c*(1+c) - (1-c)*(-k*c)) / (1+c)^2, nrow = length(x))
}

generateTrafo <- function(x0, delta, k) {

  fun <- function(x) {
    trafoFun(x, x0, delta, k)
  }

  invfun <- function(z) {
    invTrafoFun(z, x0, delta, k)
  }

  jac <- function(x) {
    trafoJac(x, x0, delta, k)
  }

  list(fun = fun, jac = jac, invfun = invfun)
}

