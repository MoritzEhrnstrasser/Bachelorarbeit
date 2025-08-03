

#### Data generation ####

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


#### save Data as single data frames for each combination ####


# Verzeichnisse vorbereiten
dir.create("sim_data", showWarnings = FALSE)

# Daten generieren und einzeln speichern
for (p in predictors) {
  for (N in N_values) {
    X <- matrix(rnorm(N * p), nrow = N, ncol = p)
    y <- rnorm(N)
    
    df <- data.frame(y = y, X)
    colnames(df) <- c("y", paste0("x", 1:p))
    
    file_name <- paste0("sim_data/data_N", N, "_p", p, ".csv")
    write.csv(df, file = file_name, row.names = FALSE)
  }
}







# Pfad zur CSV-Datenbank
data_path <- "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"  
data_files <- list.files(data_path, pattern = "*.csv", full.names = TRUE)

# Leerer Dataframe für Timing-Ergebnisse
timing_results <- data.frame(
  N = integer(),
  p = integer(),
  mean_time_us = numeric(),
  sd_time_us   = numeric()
)

# Matrix-basierte Regression mit Zeitmessung
matrix_version_timed <- function(X, y, p) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1) 
  
  t1 <- Sys.time()
  
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  
 
  y_hat <- Xmat %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - p - 1) 
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  
 
  
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = N - p - 1)
  
  t2 <- Sys.time()
  
  output_table <- cbind(beta_hat, se_beta, t_stats,p_values )
  
  as.numeric(difftime(t2, t1, units = "secs"))
}

# Benchmark-Schleife
isample <- 300

for (file in data_files) {
  df <- read.csv(file)
  
  # Y und X extrahieren
  y <- df$y
  X <- as.matrix(df[ , grep("^x", names(df))])
  
  # Metadaten aus Dateinamen
  matches <- regmatches(file, regexec("N(\\d+)_p(\\d+)", file))[[1]]
  N <- as.integer(matches[2])
  p <- as.integer(matches[3])
  
  
  # Zeitmessung über viele Wiederholungen
  times <- sapply(1:isample, function(i) matrix_version_timed(X, y, p))
  
  # Ergebnis anhängen
  timing_results <- rbind(timing_results, data.frame(
    N = N,
    p = p,
    median_time = median(times) * 1e6,
    sd_time   = sd(times) * 1e6
  ))
}

write.csv(timing_results, "timing_results_r.csv", row.names = FALSE)


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

# besser mt geom smooth






#### Julia Daten einlesen ####

julia_results <- read.csv("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/julia_whole_timed.csv")
timing_results <- read.csv("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/timing_results_r.csv")



# Beispiel: Sortieren nach den Variablen `var1` und `var2` aufsteigend
julia_results_sorted <- julia_results %>% arrange(N, p)
timing_results_sorted <- timing_results %>% arrange(N, p)


timing_results_sorted$median_time / julia_results_sorted $median_time_us

library(dplyr)
library(ggplot2)

# 2. Beide Dataframes vereinheitlichen und taggen
r_df <- timing_results %>%
  rename(median_time = median_time) %>%      
  mutate(language = "R")

julia_df <- julia_results %>%
  rename(N = N,
         p = p,
         median_time = median_time_us) %>% 
  mutate(language = "Julia")

all_df <- bind_rows(r_df, julia_df)

# 3. Gemeinsamer Plot
ggplot(all_df, aes(x = N, y = median_time, color = language, group = language)) +
  geom_line(size = 1) +
  geom_point(aes(shape = language), size = 1.5) +
  facet_wrap(~ p, scales = "free", ncol = 2,
             labeller = labeller(p = function(x) paste("p =", x))) +
  scale_color_manual(values = c(R = "#1b9e77", Julia = "#d95f02")) +
  scale_shape_manual(values = c(R = 16, Julia = 17)) +
  labs(
    title = "Vergleich der Laufzeiten von R vs. Julia\n(matrixbasierte Regression)",
    x = "Stichprobengröße N",
    y = "Median runtime (µs)",
    color = "Language",
    shape = "Languagee"
  ) +
  theme_minimal(base_size = 14)




#### geom smooth ####

ggplot(all_df, aes(x = N, y = median_time, color = language, group = language)) +
  geom_smooth(method = "loess", se = FALSE, span = 0.8, size = 1) +
  facet_wrap(~ p, scales = "free", ncol = 2,
             labeller = labeller(p = function(x) paste("p =", x))) +
  scale_color_manual(values = c(R = "#1b9e77", Julia = "#d95f02")) +
  labs(
    y = "Median runtime (µs)",
    color = "Language"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = c(0.85, 0.0001),
    legend.justification = c("right", "bottom"),
    legend.box.background = element_rect(color = "gray80", fill = "white"),
    legend.box.margin = margin(5, 5, 5, 5)
  )











#### mit Poweranalyse ####


library(tidyverse)

# ---- Einstellungen ----
data_path <- "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data"  
data_files <- list.files(data_path, pattern = "*.csv", full.names = TRUE)

isample <- 300             # Wiederholungen für Benchmark
power_iter <- 10000         # Wiederholungen für Poweranalyse
alpha <- 0.05              # Signifikanzniveau
target_N <- 3000
target_p <- 10

# Leere Tabellen vorbereiten
timing_results <- data.frame()
output_table_all <- list()
power_test_results <- NULL


# ---- Benchmarkfunktion mit Output ----
matrix_version_timed <- function(X, y, p) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1)
  
  t1 <- Sys.time()
  
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  y_hat <- Xmat %*% beta_hat
  residuals <- Y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - p - 1)
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  t_stats <- as.vector(beta_hat) / se_beta
  p_values <- 2 * pt(-abs(t_stats), df = N - p - 1)
  
  t2 <- Sys.time()
  
  time_elapsed <- as.numeric(difftime(t2, t1, units = "secs"))
  
  output_table <- data.frame(
    predictor = c("Intercept", paste0("x", 1:p)),
    beta_hat = as.vector(beta_hat),
    se = se_beta,
    t = t_stats,
    p = p_values
  )
  
  list(time = time_elapsed, output = output_table)
}


# ---- Poweranalyse für x1 ----
simulate_power <- function(N, p, beta1 = 0.0, n_iter = 10000, alpha = 0.05) {
  sig_count <- 0
  
  for (i in 1:n_iter) {
    X <- matrix(rnorm(N * p), nrow = N)
    beta <- rep(0, p)
    beta[1] <- beta1
    y <- X %*% beta + rnorm(N)
    
    Xmat <- cbind(1, X)
    XtX_inv <- solve(crossprod(Xmat))
    beta_hat <- XtX_inv %*% crossprod(Xmat, y)
    y_hat <- Xmat %*% beta_hat
    residuals <- y - y_hat
    sigma2_hat <- sum(residuals^2) / (N - p - 1)
    se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
    t_stats <- as.vector(beta_hat) / se_beta
    p_values <- 2 * pt(-abs(t_stats), df = N - p - 1)
    
    if (p_values[2] < alpha) sig_count <- sig_count + 1  # x1 nach Intercept
  }
  
  sig_count / n_iter
}


# ---- Hauptschleife ----
for (file in data_files) {
  df <- read.csv(file)
  y <- df$y
  X <- as.matrix(df[ , grep("^x", names(df))])
  
  matches <- regmatches(file, regexec("N(\\d+)_p(\\d+)", file))[[1]]
  N <- as.integer(matches[2])
  p <- as.integer(matches[3])
  
  cat("Verarbeite Datei: N =", N, "p =", p, "\n")
  
  times <- numeric(isample)
  outputs <- list()
  
  for (i in 1:isample) {
    result <- matrix_version_timed(X, y, p)
    times[i] <- result$time
    
    # Output-Tabelle nur für die Zielkombination speichern
    if (N == target_N && p == target_p && i == 1) {
      outputs[[i]] <- result$output
    }
  }
  
  # Mittelwert-Output speichern
  if (N == target_N && p == target_p) {
    output_table_all <- result$output
    
    # Poweranalyse durchführen
    power_val <- simulate_power(N = N, p = p, beta1 = 0,
                                n_iter = power_iter, alpha = alpha)
    power_test_results <- data.frame(N = N, p = p, power = power_val, alpha = alpha)
  }
  
  # Timing speichern
  timing_results <- rbind(timing_results, data.frame(
    N = N,
    p = p,
    median_time_us = median(times) * 1e6,
    sd_time_us     = sd(times) * 1e6
  ))
}


# ---- Export ----
write.csv(timing_results, "timing_results_r.csv", row.names = FALSE)
write.csv(output_table_all, "output_table_N3000_p10.csv", row.names = FALSE)

if (!is.null(power_test_results)) {
  write.csv(power_test_results, "power_analysis_N3000_p10.csv", row.names = FALSE)
}

cat("Fertig! Ergebnisse gespeichert:\n")
cat("• timing_results_r.csv\n")
cat("• output_table_N3000_p10.csv\n")
cat("• power_analysis_N3000_p10.csv\n")




df <- read.csv("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data/data_N2000_p5.csv")

fit_alm <- lm(y~ x1 + x2 + x3 + x4 + x5, data = df)

summary(fit_alm)




# comparing theoretical and empirical SE

# Simulationsparameter
n_sim <- 10000  # Anzahl der Simulationen
N <- 1000        # Stichprobengröße
p <- 3          # Anzahl Prädiktoren
beta <- c(1, 0.5, -0.3, 0.2)  # Wahrer Koeffizientenvektor (inkl. Intercept)
sigma <- 1

# Ergebnisse speichern
estimates <- matrix(NA, nrow = n_sim, ncol = length(beta))
se_theoretical <- matrix(NA, nrow = n_sim, ncol = length(beta))

for (i in 1:n_sim) {
  # Simuliere Daten
  X <- matrix(rnorm(N * p), ncol = p)
  Xmat <- cbind(1, X)
  y <- Xmat %*% beta + rnorm(N, sd = sigma)
  
  # OLS-Schätzung via Matrix-Algebra
  XtX_inv <- solve(crossprod(Xmat))
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  y_hat <- Xmat %*% beta_hat
  residuals <- y - y_hat
  sigma2_hat <- sum(residuals^2) / (N - p - 1)
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))
  
  # Speichern
  estimates[i, ] <- as.vector(beta_hat)
  se_theoretical[i, ] <- se_beta
}

# Empirischer Standardfehler = SD der Schätzungen
se_empirical <- apply(estimates, 2, sd)
se_theoretical_mean <- colMeans(se_theoretical)

# Vergleich in DataFrame
comparison_df <- data.frame(
  Term = paste0("β", 0:p),
  Empirical_SE = round(se_empirical, 4),
  Mean_Theoretical_SE = round(se_theoretical_mean, 4)
)

print(comparison_df)




