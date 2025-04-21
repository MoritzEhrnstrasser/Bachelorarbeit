

library(microbenchmark)


# ---- ALM Output replizieren mit matrixschreibweise ----

isample = 10000
N = 2500
B0 = 0
B1 = 0
sd = 1

# Generate data
generate_data<-function(N, B0, B1, sd) {
  x <- rnorm(N, mean=0, sd = 1)
  e <- rnorm(N , mean = 0, sd = 1)
  e <- e*(sqrt(1)/sd(e))
  y <- B0 + B1 * x + e
  
  fit <- lm(y ~ x)
  sum_fit <-  summary(fit)$coefficients
  
  p_value <- sum_fit["x","Pr(>|t|)"]
  return(p_value)
}

p_values <- sapply(1:isample, function(i) generate_data(N, B0, B1, sd))

mean(p_values < 0.05) # 0.05 is the significance level





# ---- ALM Output replizieren mit Matrixschreibweise ----

isample = 10000
N = 2500
B0 = 0
B1 = 0
sd = 1

# Generate data
generate_data <- function(N, B0, B1, sd) {
  x <- rnorm(N, mean = 0, sd = 1)
  e <- rnorm(N, mean = 0, sd = 1)
  e <- e * (sqrt(1) / sd(e))  # Normalize error
  y <- B0 + B1 * x + e
  
  # Matrix form
  X <- cbind(1, x)              # Design matrix with intercept
  Y <- matrix(y, ncol = 1)      # Make Y a column vector
  
  # Beta estimate: (X'X)^(-1) X'Y
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% Y  # 2x1 matrix
  
  # Residuals and variance estimate
  y_hat <- X %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  
  # Standard errors: sqrt(diag(sigma² * (X'X)^(-1)))
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  
  # t-statistics
  t_stats <- as.vector(beta_hat) / se_beta
  
  # p-values (2-sided)
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  output_table <- cbind(beta_hat, se_beta, t_stats,p_values )
  
  
  fit <- lm( y ~ x)
  summary(fit)
  
  return(p_values[2])  # Return p-value of x (not intercept)
}

# Monte Carlo simulation
p_values <- sapply(1:isample, function(i) generate_data(N, B0, B1, sd))

# Proportion of significant p-values
mean(p_values < 0.05)



# ---- Benchmarking speed ----
#### daten erstellen ####
N <- 250
B0 = 0
B1 = 0
x <- rnorm(N)
e <- rnorm(N)
y <- B0 + B1 * x + e


# lm()-Variante
lm_version <- function() {
  fit <- lm(y ~ x)
  summary(fit)
}


# Matrix-Version
matrix_version <- function() {
  X <- cbind(1, x)
  Y <- matrix(y, ncol = 1)      # Make Y a column vector
  
  # Beta estimate: (X'X)^(-1) X'Y
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% Y  # 2x1 matrix
  
  # Residuals and variance estimate
  y_hat <- X %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  
  # Standard errors: sqrt(diag(sigma² * (X'X)^(-1)))
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  
  # t-statistics
  t_stats <- as.vector(beta_hat) / se_beta
  
  # p-values (2-sided)
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  output_table <- cbind(beta_hat, se_beta, t_stats,p_values )
  output_table
}

microbenchmark(
  lm = lm_version(),
  matrix = matrix_version(),
  times = 100L  # 100 Wiederholungen
)





# ---- Benchmarking speed depending on  sample size ----
install.packages(c("microbenchmark", "ggplot2"))   # falls noch nicht installiert
library(microbenchmark)
library(ggplot2)

# Funktionen definieren
lm_version <- function(x, y) {
  fit <- lm(y ~ x)
  summary(fit)
}

matrix_version <- function(x, y) {
  N <- length(y)
  X <- cbind(1, x)
  Y <- matrix(y, ncol = 1)      # Make Y a column vector
  
  # Beta estimate: (X'X)^(-1) X'Y
  XtX_inv <- solve(t(X) %*% X)
  beta_hat <- XtX_inv %*% t(X) %*% Y  # 2x1 matrix
  
  # Residuals and variance estimate
  y_hat <- X %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  
  # Standard errors: sqrt(diag(sigma² * (X'X)^(-1)))
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  
  # t-statistics
  t_stats <- as.vector(beta_hat) / se_beta
  
  # p-values (2-sided)
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  output_table <- cbind(beta_hat, se_beta, t_stats,p_values )
  output_table
}

# Vektor von Stichproben-Größen
N_values <- seq(250, 5000, by = 100)

# Ergebnis-Speicher
results <- data.frame(
  N = integer(),
  method = character(),
  time_us = numeric(),
  stringsAsFactors = FALSE
)

# ---- Benchmark-Schleife ----
set.seed(42)
for (N in N_values) {
  # Einmal Daten erzeugen
  x <- rnorm(N)
  y <- 0 + 0 * x + rnorm(N)
  
  bm <- microbenchmark(
    lm     = lm_version(x, y),
    matrix = matrix_version(x, y),
    times  = 50L
  )
  
  # Median-Zeit in Mikrosekunden
  mb <- summary(bm)[, c("expr", "median")]
  colnames(mb) <- c("method", "time_us")
  mb$N <- N
  
  # Anfügen
  results <- rbind(results, mb[, c("N", "method", "time_us")])
}

# ---- Plot ----
ggplot(results, aes(x = N, y = time_us, color = method)) +
  geom_line(size = 1) +
  geom_point(size = 1.5) +
  scale_y_continuous(name = "Median-Laufzeit (µs)") +
  scale_x_continuous(name = "Stichprobengröße N") +
  ggtitle("Laufzeit von lm() vs. Matrix-Regression") +
  theme_minimal(base_size = 14)


