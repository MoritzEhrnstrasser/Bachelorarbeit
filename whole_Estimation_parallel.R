


#### sequenziell vs parallel ####



# ==== R-Skript: Seq. vs. Parallel mit future.apply ====

library(future.apply)

# 1. Setup
data_path   <- "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"
data_files  <- list.files(data_path, pattern = "\\.csv$", full.names = TRUE)
isample     <- 3000
n_cores     <- 4

# 2. Regressionsfunktion
matrix_version_timed <- function(X, y, p) {
  t1 <- Sys.time()
  Xmat <- cbind(1, X)
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  y_hat <- Xmat %*% beta_hat
  residuals <- y - y_hat
  sigma2_hat <- sum(residuals^2) / (length(y) - p - 1)
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = length(y) - p - 1)
  t2 <- Sys.time()
  as.numeric(difftime(t2, t1, units = "secs"))
}

# 3. Datei-Benchmark-Funktion
benchmark_single_file <- function(file) {
  df <- read.csv(file)
  y  <- df$y
  X  <- as.matrix(df[, grep("^x", names(df))])
  # N und p aus Dateiname
  caps <- regmatches(file, regexec("N(\\d+)_p(\\d+)", file))[[1]]
  N    <- as.integer(caps[2])
  p    <- as.integer(caps[3])
  # Messung
  times <- replicate(isample, matrix_version_timed(X, y, p))
  data.frame(N = N,
             p = p,
             median_time_us = median(times) * 1e6,
             sd_time_us     = sd(times)     * 1e6)
}

# 4.a Sequentiell
start_seq <- Sys.time()
results_seq <- do.call(rbind, lapply(data_files, benchmark_single_file))
end_seq   <- Sys.time()
t_seq     <- as.numeric(difftime(end_seq, start_seq, units = "secs"))

# 4.b Parallel
plan(multisession, workers = n_cores)
start_par <- Sys.time()
results_par <- future_sapply(data_files, benchmark_single_file, simplify = FALSE)
results_par <- do.call(rbind, results_par)
end_par   <- Sys.time()
t_par     <- as.numeric(difftime(end_par, start_par, units = "secs"))

# 5. Speedup & Effizienz
speedup    <- t_seq / t_par
efficiency <- speedup / n_cores


cat("R Benchmark (n_cores =", n_cores, ")\n")
cat(" Sequentiell:", round(t_seq, 2), "s\n")
cat(" Parallel:   ", round(t_par, 2), "s\n")
cat(" Speedup:    ", round(speedup, 2), "\n")
cat(" Effizienz:  ", round(efficiency * 100, 1), "%\n")


# 6. Export optional
write.csv(results_seq, "/tmp/timing_seq_r.csv", row.names = FALSE)
write.csv(results_par, "/tmp/timing_par_r.csv", row.names = FALSE)


