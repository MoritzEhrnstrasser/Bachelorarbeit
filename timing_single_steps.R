library(dplyr)
library(tidyverse)

matrix_steps_timed <- function(X, y,p) {
  N <- length(y)
  Xmat <- cbind(1, X)
  Y <- matrix(y, ncol = 1)
  
  t1 <- Sys.time()
  XtX_inv <- solve(crossprod(Xmat))               # Inversion
  beta_hat <- XtX_inv %*% crossprod(Xmat, y)
  t2 <- Sys.time()
  
  y_hat <- Xmat %*% beta_hat                      # Vorhersage
  residuals <- Y - y_hat
  t3 <- Sys.time()
  
  sigma2_hat <- sum(residuals^2) / (N - p - 1)
  
  se_beta <- sqrt(diag(sigma2_hat * XtX_inv))     # Standardfehler
  t4 <- Sys.time()
  
  t_stats <- as.vector(beta_hat) / se_beta        # p-Werte
  p_values <- 2 * pt(-abs(t_stats), df = N - 2)
  t5 <- Sys.time()
  
  c(
    inversion = as.numeric(difftime(t2, t1, units = "secs")),
    prediction = as.numeric(difftime(t3, t2, units = "secs")),
    standard_error = as.numeric(difftime(t4, t3, units = "secs")),
    p_values = as.numeric(difftime(t5, t4, units = "secs"))
  )
}

# Benchmark-Durchlauf
timing_steps_df <- data.frame(N = integer(), p = integer(),
                              inversion = numeric(), prediction = numeric(),
                              standard_error = numeric(), p_values = numeric())

# Lade z.B. CSV-Dateien
data_files <- list.files("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_data", full.names = TRUE)

for (file in data_files) {
  df <- read.csv(file)
  y <- df$y
  X <- as.matrix(df[ , grep("^x", names(df))])
  
  matches <- regmatches(file, regexec("N(\\d+)_p(\\d+)", file))[[1]]
  N <- as.integer(matches[2])
  p <- as.integer(matches[3])
  
  step_times <- replicate(500, matrix_steps_timed(X, y,p))
  step_medians <- apply(step_times, 1, median)  
  step_medians <- step_medians * 1e6
  
  timing_steps_df <- rbind(timing_steps_df,
                           data.frame(N = N, p = p, t(step_medians)))
}

write.csv(timing_steps_df, "timing_steps_r.csv", row.names = FALSE)






####  fertiger Plot ####

library(dplyr)
library(tidyr)
library(ggplot2)

# Einlesen und Label
r <- read.csv("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/timing_steps_r.csv") %>% 
  mutate(lang = "R")
j <- read.csv("/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/timing_steps_julia.csv") %>% 
  mutate(lang = "Julia")

# Zusammenführen
all_steps <- bind_rows(r, j)

# Long Format + Umbenennung
long_steps <- all_steps %>%
  pivot_longer(cols = c(inversion, p_values, standard_error, prediction),
               names_to = "step", values_to = "time_us") %>%
  mutate(step = recode(step,
                       inversion = "inversion and beta estimation",
                       p_values = "p values",
                       standard_error = "standard error"))

# Stacked Area Plot
ggplot(long_steps, aes(x = N, y = time_us, fill = step)) +
  geom_area(position = "stack", alpha = 0.8) +
  facet_grid(lang ~ p, labeller = label_both, scales = "free_y") +
  labs(
    x = "Sample size N",
    y = "Runtime (µs)",
    fill = "Estimation step"
  ) +
  theme_minimal(base_size = 16) +  # Größere Basisschrift
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 9),
    strip.text = element_text(size = 16),      # Facet-Titel
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    plot.title = element_text(size = 20, face = "bold")
  )
+
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "bottom")


# vertikal angeordnet

long_steps$lang <- factor(long_steps$lang, levels = c("R", "Julia"))

ggplot(long_steps, aes(x = N, y = time_us, fill = step)) +
  geom_area(position = "stack", alpha = 0.9) +
  facet_grid(p ~ lang, labeller = label_both, scales = "free_y") +
  
  labs(
    x = "Sample size N",
    y = "Runtime (µs)",
    fill = "Estimation step"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 13),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines")
  )



