library(Matrix)

corr_func <- function(x, length_scale) {
  d_over_l <- abs(outer(x,x,"-"))/length_scale
  
  idx <- which(d_over_l<1, arr.ind=TRUE)
  d_over_l <- d_over_l[idx]
  data <- (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)
  
  sparseMatrix(i=idx[,1], j=idx[,2], x=data)
}

x <- seq(0,3,1)
corr_mat <- corr_func(x, length_scale=5.0)

sigmas <- 10*runif(length(x))

# first way
D <- diag(sigmas,nrow=length(x))
S1 <- D %*% corr_mat %*% t(D) # covaraince matrix

# second way
S2 <- (corr_mat * sigmas) * rep(sigmas, each=nrow(corr_mat))

# third way
S3 <- outer(sigmas,sigmas) * corr_mat

####################################3

x <- seq(1,3,1)

A <- matrix(seq(1,9), nrow=3)

x_xt <- outer(x, x)

x_xt*A
