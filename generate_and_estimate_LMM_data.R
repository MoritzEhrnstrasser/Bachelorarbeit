
# 1. Benötigte Pakete
library(Matrix)    # für bdiag
library(MASS)      # für mvrnorm

# 2. Simulationsfunktion
simulate_lmm_data <- function(model_type = 1,
                              p,
                              n_per_cluster,
                              n_clusters,
                              beta = rep(0, p),
                              sigma_u = 1.5,
                              sigma_eps = 1,
                              seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  N <- n_per_cluster * n_clusters
  # Design X: Intercept + p Prädiktoren
  X_pred <- matrix(rnorm(N * p), nrow = N, ncol = p)
  X <- cbind(1, X_pred)
  colnames(X) <- c("Intercept", paste0("x", 1:p))
  
  # Cluster-IDs
  cluster_id <- rep(1:n_clusters, each = n_per_cluster)
  
  # Zufallseffekte generieren
  if (model_type == 1) {
    # Modell 1: nur random intercept
    Z <- Matrix::bdiag(lapply(1:n_clusters, function(i) matrix(1, n_per_cluster, 1)))
    u <- rnorm(n_clusters, 0, sqrt(sigma_u))
    
  } else {
    # Modell 2: random intercept + random slopes für alle Prädiktoren
    Zlist <- lapply(1:n_clusters, function(i) {
      rows <- ((i-1)*n_per_cluster + 1):(i*n_per_cluster)
      as.matrix(X[rows, ])  # (n_per_cluster) × (p+1)
    })
    Z <- Matrix::bdiag(Zlist)
    # multivariate Zufallseffekte: Intercept + p Slopes
    u <- as.vector(MASS::mvrnorm(n_clusters, mu = rep(0, p+1),
                                 Sigma = diag(sigma_u, p+1)))
  }
  
  # Residualfehler 
  e <- rnorm(N, 0, sqrt(sigma_eps))
  
  # Linearkombination
  y_fixed  <- X %*% c(0, beta)      # beta ohne zusätzlichen intercept-Koeff
  y_random <- as.vector(Z %*% u)
  y        <- y_fixed + y_random + e
  
  # Rückgabe als data.frame
  df <- data.frame(
    y = y,
    cluster_id = factor(cluster_id),
    X_pred
  )
  names(df)[3:(3+p-1)] <- paste0("x", 1:p)
  return(df)
}

# 3. Parameter-Bereiche
ps         <- c(3, 6)
clusters   <- seq(200, 400, by = 50)
per_clust  <- seq(200, 400, by = 50)
out_dir    <- "sim_lmm_data"
dir.create(out_dir, showWarnings = FALSE)

# 4. Schleifen und Speichern
for (model in 1:2) {
  # model==1: Random Intercept; model==2: Random Intercept + Slopes
  model_label <- ifelse(model == 1, "ri", "rs")
  
  for (p in ps) {
    for (n_c in clusters) {
      for (n_pc in per_clust) {
        df <- simulate_lmm_data(
          model_type      = model,
          p               = p,
          n_per_cluster   = n_pc,
          n_clusters      = n_c,
          beta            = rep(0.3, p),
          sigma_u         = 1.5,
          sigma_eps       = 1,
          cor_u = 0
        )
        
        file_name <- sprintf("%s/lmm_%s_p%d_nc%d_npc%d.csv",
                             out_dir, model_label, p, n_c, n_pc)
        write.csv(df, file = file_name, row.names = FALSE)
      }
    }
  }
}

cat("✅ Fertig: Daten in", out_dir, "gespeichert.\n")






#### bereinigter Code ohne Doppelungen ####


library(lme4)
library(dplyr)
library(ggplot2)

# Pfad zu den simulierten LMM-Daten
data_path <- "/Users/moritz/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/sim_lmm_data" 
data_files <- list.files(data_path, pattern = "\\.csv$", full.names = TRUE)

# Benchmark-Einstellungen
isample <- 10    # Anzahl der Wiederholungen pro Datei

# Ergebnis-Container
results <- data.frame(
  model_type = character(),
  p          = integer(),
  n_clusters = integer(),
  n_per_cl   = integer(),
  median_us  = numeric(),
  sd_us      = numeric(),
  stringsAsFactors = FALSE
)

# Hilfsfunktion: Zeitmessung mit Sys.time()
time_lmer <- function(formula, data, reps = isample) {
  times <- numeric(reps)
  for (i in seq_len(reps)) {
    t_start <- Sys.time()
    lmer(formula, data = data, REML = TRUE)
    t_end <- Sys.time()
    times[i] <- as.numeric(difftime(t_end, t_start, units = "secs"))
  }
  times * 1e6  # in Mikrosekunden
}

# Loop über alle Dateien
for (file in data_files) {
  df <- read.csv(file)
  
  # Metadaten aus Dateiname extrahieren
  parts <- strcapture(
    "lmm_(ri|rs)_p([0-9]+)_nc([0-9]+)_npc([0-9]+)\\.csv",
    basename(file),
    proto = list(model=character(), p=integer(), nc=integer(), npc=integer())
  )
  model_label <- parts$model       # "ri" oder "rs"
  p_val       <- parts$p
  nc          <- parts$nc
  npc         <- parts$npc
  
  # Formeln bauen
  fixed_part <- paste0("x", 1:p_val, collapse = " + ")
  formula_ri <- as.formula(paste("y ~", fixed_part, "+ (1 | cluster_id)"))
  formula_rs <- as.formula(
    paste0("y ~ ", fixed_part, " + (1 + ", fixed_part, " | cluster_id)")
  )
  
  # Passendes Modell auswählen
  if (model_label == "ri") {
    times <- time_lmer(formula_ri, df, isample)
  } else {
    times <- time_lmer(formula_rs, df, isample)
  }
  
  # Ergebnis speichern
  summary_df <- data.frame(
    model_type = model_label,
    p          = p_val,
    n_clusters = nc,
    n_per_cl   = npc,
    median_us  = median(times),
    sd_us      = sd(times)
  )
  
  results <- bind_rows(results, summary_df)
}

# Ergebnisse speichern
write.csv(results, "benchmark_lmm_lme4_sys_time_clean.csv", row.names = FALSE)

r_lmm <- read.csv("~/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/benchmark_lmm_lme4_sys_time_clean.csv")

# Beispielplot 
r_lmm %>%
  filter(p == 6) %>%
  ggplot(aes(x = n_clusters * n_per_cl, y = median_us/1e6,
             color = model_type)) +
  geom_line() +
  labs(title = "predictor size = 6",
       x = "total sample size = n cluster x n per cluster",
       y = "median runtime (s)",
       color = "model") +
  theme_minimal()


# Beispielplot 
r_lmm %>%
  filter(p == 3) %>%
  ggplot(aes(x = n_clusters * n_per_cl, y = median_us/1e6,
             color = model_type)) +
  geom_line() +
  labs(title = "predictor size = 6",
       x = "total sample size = n cluster x n per cluster",
       y = "median runtime (s)",
       color = "model") +
  theme_minimal()


# Julia Daten einlesen
j_lmm <- read.csv("~/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/benchmark_lmm_mixedmodels.csv")

j_lmm %>%
  filter(p == 3) %>%
  ggplot(aes(x = n_clusters * n_per_cl, y = median_us/1e6,
             color = model_type)) +
  geom_line() +
  labs(title = "predictor size = 3",
       x = "total sample size = n cluster x n per cluster",
       y = "median runtime (s)",
       color = "model") +
  theme_minimal()

j_lmm %>%
  filter(p == 6) %>%
  ggplot(aes(x = n_clusters * n_per_cl, y = median_us/1e6,
             color = model_type)) +
  geom_line() +
  labs(title = "predictor size = 3",
       x = "total sample size = n cluster x n per cluster",
       y = "median runtime (s)",
       color = "model") +
  theme_minimal()





#### plot als grid ####
library(dplyr)
library(ggplot2)
library(readr)

# R-Daten laden und Sprache labeln
r_lmm <- read_csv("~/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/benchmark_lmm_lme4_sys_time_clean.csv") %>%
  mutate(lang = "R")

# Julia-Daten laden und Sprache labeln
j_lmm <- read_csv("~/Documents/Uni Kram/R/Projekte/Bachelorarbeit/Bachelorarbeit R/results/benchmark_lmm_mixedmodels.csv") %>%
  mutate(lang = "Julia")

# Kombinieren
all_lmm <- bind_rows(r_lmm, j_lmm)


all_lmm$lang <- factor(all_lmm$lang, levels = c("R", "Julia"))
# Plot mit facet_grid für Sprache und Prädiktoranzahl
ggplot(all_lmm, aes(x = n_clusters * n_per_cl, y = median_us / 1e6, color = model_type)) +
  geom_smooth(size = 0.9) +
  facet_grid(lang ~ p, labeller = label_both, scales = "free_y") +
  labs(
    
    x = expression("Total sample size (" * n[clusters] * " × " * n[per~cluster] * ")")
    ,
    y = "Median runtime (s)",
    color = "Model type"
  ) +
  
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = "bottom",
    legend.key.size = unit(0.8, "lines")
  )

citation("afex")
