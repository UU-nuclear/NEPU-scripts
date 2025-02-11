cov_func <- function(d, length_scale) {
	d_over_l <- abs(d)/length_scale
	ifelse(d_over_l < 1, ( (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)), 0)
	# OBS: The sigma hyper-parameter is not included here but outside the function
}

matern_1_2 <- function(d, length_scale) {
		exp(-abs(d)/length_scale)
}

matern_3_2 <- function(d, length_scale) {
	(1+sqrt(3)*abs(d)/length_scale)*exp(-abs(x1-x2)/length_scale)
}

sqr_exp <- function(d, length_scale) {
	exp(-0.5*d^2/length_scale^2)
}

x <-  seq(from=-5, to=5, by=0.001)
y_sparse <- cov_func(x, length_scale=3.121378)
y_mat_1_2 <- matern_1_2(x, length_scale=1/log(2))
y_prod <- y_sparse*y_mat_1_2

plot(x,y_mat_1_2, type='l')
lines(x[y_sparse>0],y_sparse[y_sparse>0], col='red')
lines(x[y_prod>0],y_prod[y_prod>0], col='green')

#=============================

y_sparse <- cov_func(x, length_scale=1)
y_mat_1_2 <- matern_1_2(x, length_scale=1/3.121378/log(2))
y_prod <- y_sparse*y_mat_1_2

plot(x,y_mat_1_2, type='l')
lines(x[y_sparse>0],y_sparse[y_sparse>0], col='red')
lines(x[y_prod>0],y_prod[y_prod>0], col='green')