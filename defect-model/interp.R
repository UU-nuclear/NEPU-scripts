library(Matrix)

# Function to create a sparse matrix for linear interpolation
interpolation_matrix <- function(x, xi) {

 # first strip of the values in xi that are out of the range of x
 xi_inside <- xi[xi>=min(x) & xi<=max(x)]

 # count the number of striped values to the left of min(x)
 n_below <- length(xi[xi<min(x)])

 dx <- diff(x)
 if (any(dx == 0)) {
   stop("creation of interpolation_matrix: x values must be distinct")
 }

 # Find the indices for xi in x
 ix <- findInterval(xi_inside, x, rightmost.closed=TRUE)

 # Calculate weights for interpolation
 wx <- (xi_inside - x[ix]) / dx[ix]
 #wx <- pmax(pmin(wx, 1), 0)  # Ensure weights are between 0 and 1

 # Create the sparse interpolation matrix
 matrix_rows <- length(xi)
 matrix_cols <- length(x)

 # each row has two entries, shift the rows to the left by the number of points out of range
 rows <- rep(seq_len(length(xi_inside)), each = 2) + n_below
 # the coloumns correspond to the indicies of x which are interpolated between
 cols <- as.vector(rbind(ix, ix+1))
 # the values is the partial derivative of the interpolated function value
 # with respect to the y-value at the energy grid x
 values <- as.vector(rbind(1-wx, wx))

 interpolation_matrix <- sparseMatrix(
   i = rows,
   j = cols,
   x = values,
   dims = c(matrix_rows, matrix_cols)
 )

 return(interpolation_matrix)

}


# Example usage
#x <- c(0,sort(runif(3,1,5)),6)
x <- 1:5
y <- runif(length(x),0,1)
# xi <- sort(runif(3,1,5))  # Values to interpolate
# xi <- runif(30,0,5)  # Values to interpolate
xi <- seq(from=0, to=6, by=0.5)  # Values to interpolate

sparse_interpolation_matrix <- interpolation_matrix(x, xi)
print(as.matrix(sparse_interpolation_matrix))

yi <- sparse_interpolation_matrix %*% y

plot(x, y, type="l")
points(xi,yi)


plot(xi, yi, type="p")
lines(x,y)

####################

sparse_matrix_list <- list(
  sparseMatrix(i = c(1, 2, 3), j = c(2, 1, 3), x = c(4, 2, 6)),  # Sparse matrix 1
  sparseMatrix(i = c(1, 2), j = c(2, 3), x = c(5, 7))            # Sparse matrix 2
)

# Combine matrices along the diagonal
result_matrix <- bdiag(sparse_matrix_list)

print(result_matrix)

######################

# Sample function
my_function <- function(arg1, arg2, arg3) {
  print(paste("Arguments:", arg1, arg2, arg3))
}

# Sample lists
my_list <- c(1, 2, 3)
arg2_values <- c("SecondArg1", "SecondArg2", "SecondArg3")
arg3_values <- c("ThirdArg1", "ThirdArg2", "ThirdArg3")

result <- Map(my_function, my_list, arg2_values, arg3_values)

result2 <- mapply(my_function, my_list, arg2_values, arg3_values, SIMPLIFY=FALSE)

######################

# Sample function
my_function <- function(arg1, arg2, arg3) {
  print(paste("Arguments:", arg1, arg2, arg3))
}

# Sample lists
my_list <- c(1, 3)
arg2_values <- c("SecondArg1", "SecondArg3")
arg3_values <- c("ThirdArg1", "ThirdArg3")

result2 <- mapply(my_function, my_list, arg2_values, arg3_values, SIMPLIFY=FALSE, MoreArgs)