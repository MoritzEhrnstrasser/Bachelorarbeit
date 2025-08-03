#### #### 

library(ggplot2)

# Parameter
predictors <- c(1, 5, 10, 20, 30)
N_values <- seq(1000, 10000, by = 250)
isample <- 500

# 1. Datensätze für jede Kombination von N und p vorbereiten und speichern
dataset_list <- list()

for (p in predictors) {
  for (N in N_values) {
    
    X <- matrix(rnorm(N * p), nrow = N, ncol = p)
    y <- rnorm(N)
    
    dataset_list[[paste0("N", N, "_p", p)]] <- list(X = X, y = y, N = N, p = p)
  }
}

# 2. Matrix-basierte Regressionsfunktion mit Zeitmessung
matrix_version_timed <- function(X, y) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1) 
  
  t3 <- Sys.time()
  t1 <- Sys.time()
  
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  y_hat <- Xmat %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  t2 <- Sys.time()
  t4 <- Sys.time()
  as.numeric(difftime(t4, t3, units = "secs"))  # Rückgabe: Zeit in Sekunden
}

# 3. Für jede Kombination Zeitmessung durchführen
timing_results <- data.frame(N = integer(),
                             p = integer(),
                             mean_time_us = numeric(),
                             sd_time_us   = numeric())

for (key in names(dataset_list)) {
  dat <- dataset_list[[key]]
  
  times <- sapply(1:isample, function(i) matrix_version_timed(dat$X, dat$y))
  
  timing_results <- rbind(timing_results, data.frame(
    N = dat$N,
    p = dat$p,
    mean_time = mean(times) * 1e6,  # µs
    sd_time   = sd(times) * 1e6
  ))
}



### Plot ####

ggplot(timing_results, aes(x = N, y = mean_time)) +
  geom_line(color = "blue") +
  geom_point(size = 0.8, color = "blue") +
  facet_wrap(~ p, scales = "free", ncol = 2,
             labeller = labeller(p = function(x) paste("p =", x))) +
  labs(
    title = "Matrixbasierte Regression: Laufzeit in Abhängigkeit von N und p",
    x = "Stichprobengröße N",
    y = "Mittlere Laufzeit (µs)"
  ) +
  theme_minimal(base_size = 14)


#### measuring Sys.time overhead ####



# Ohne Zeitmessung
matrix_version_plain <- function(X, y) {
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
  
  invisible(NULL)
}

# Mit Zeitmessung (wie bisher)
matrix_version_timed <- function(X, y) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1)
  
  t1 <- Sys.time()
  
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  y_hat <- Xmat %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - 2)
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  
  t2 <- Sys.time()
  as.numeric(difftime(t2, t1, units = "secs"))
}



# Beispiel-Daten

N <- 5000
p <- 10
X <- matrix(rnorm(N * p), nrow = N, ncol = p)
y <- rnorm(N)

# Wiederholungen
reps <- 1000

# Zeitmessung für beide Varianten
time_plain <- system.time({
  for (i in 1:reps) matrix_version_plain(X, y)
})

time_timed <- system.time({
  for (i in 1:reps) matrix_version_timed(X, y)
})


