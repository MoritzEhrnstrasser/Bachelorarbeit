# ---- microbenchmark vs sys.time ----


library(microbenchmark)

# estimate regression
regression_matrix <- function(X, y) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1) 
  
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  
  y_hat <- Xmat %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  output_table <- cbind(beta_hat, se_beta, t_stats, p_values)
  return(output_table)
}

# parameter
isample <- 1000
N <- 2500
p <- 10

# generate Data
X <- matrix(rnorm(N * p), nrow = N)
y <- rnorm(N)

# measuring Sys.time 
times_sys <- sapply(1:isample, function(i) {
  t1 <- Sys.time()
  regression_matrix(X, y)
  t2 <- Sys.time()
  as.numeric(difftime(t2, t1, units = "secs"))
})


median_time <- median(times_sys)
sd_time <- sd(times_sys)

cat("  Median:", format(mean_time * 1e6, digits = 4), "µs\n")
cat("  SD:       ", format(sd_time  * 1e6, digits = 4), "µs\n")


# 4) measuring Sys.time
mb <- microbenchmark(
  regression_matrix(X, y),
  times = isample,
  unit  = "us"    # in microseconds
)

summary(mb)


