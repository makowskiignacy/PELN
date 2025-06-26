library(tidyverse)
library(janitor)
library(quantreg)

# Load the datasets, specifying comma as decimal mark
ref_df <- read_csv('data/PELN - Van Wagner całkowity pomiar (martwe drewno).csv', locale = locale(decimal_mark = ","))
wz_df <- read_csv('data/PELN - Van Wagner WZ (martwe drewno).csv', locale = locale(decimal_mark = ","))
ns_df <- read_csv('data/PELN - Van Wagner NS (martwe drewno).csv', locale = locale(decimal_mark = ","))

# Clean column names to be more script-friendly
ref_df <- ref_df %>% clean_names()
wz_df <- wz_df %>% clean_names()
ns_df <- ns_df %>% clean_names()

# --- Constants ---
AREA_HA <- 0.5
L_PER_TRANSECT <- 10 # Assuming each transect line is 10m long

# --- Total Reference Volume Calculation (per Hectare) ---
total_ref_vol_m3 <- ref_df %>%
  mutate(
    diameter_m = srednica_w_polowie_dlugosci_cm / 100,
    length_m = dlugosc_m
  ) %>%
  mutate(volume_m3 = pi * (diameter_m / 2)^2 * length_m) %>%
  summarise(total_volume = sum(volume_m3, na.rm = TRUE)) %>%
  pull(total_volume)

# Reference volume per hectare
ref_vol_m3_ha <- total_ref_vol_m3 / AREA_HA

# --- NS Volume Estimation (per Square) ---
# The Van Wagner formula V(m³/ha) = (π² * Σdᵢ²[cm²]) / (8 * L[m])
ns_estimates <- ns_df %>%
  mutate(d2 = srednica_w_przecieciu^2) %>%
  group_by(nr_kwadratu) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT))

# --- WZ Volume Estimation (per Line) ---
wz_estimates <- wz_df %>%
  mutate(d2 = srednica_w_przecieciu^2) %>%
  group_by(nr_linii) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT))

# --- Data Preparation for Plotting ---
plot_df <- bind_rows(
  ns_estimates %>% select(volume_m3_ha) %>% mutate(method = "NS per Square"),
  wz_estimates %>% select(volume_m3_ha) %>% mutate(method = "WZ per Line")
)

# --- Plotting ---
# Define a more elegant color palette
elegant_palette <- c("NS per Square" = "#004d00", "WZ per Line" = "#697a87") # Dark green and slate grey

violin_plot <- ggplot(plot_df, aes(x = method, y = volume_m3_ha)) +
  stat_ydensity(geom = "violin", aes(fill = method), trim = FALSE, alpha = 0.6, bounds = c(0, Inf)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  # Add a point marker for the reference value on each violin plot
  geom_point(data = . %>% distinct(method), aes(y = ref_vol_m3_ha, shape = "Reference"), size = 5, color = "black") +
  scale_fill_manual(values = elegant_palette, name = "Method") +
  scale_shape_manual(name = "", values = c("Reference" = 18), labels = c(sprintf("Reference: %.2f m³/ha", ref_vol_m3_ha))) +
  labs(
    title = "Distribution of Volume Estimates vs. Reference Value",
    subtitle = "Estimates calculated per NS Square and WZ Line",
    x = "Method",
    y = "Volume (m³/ha)",
    fill = "Method"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("volume_distribution_vs_reference.png", violin_plot, width = 8, height = 7)

# --- Statistical Significance ---
# One-sample t-test to compare the mean of estimates to the reference value
ns_ttest <- t.test(ns_estimates$volume_m3_ha, mu = ref_vol_m3_ha)
wz_ttest <- t.test(wz_estimates$volume_m3_ha, mu = ref_vol_m3_ha)

cat("--- Statistical Analysis ---
")
cat(sprintf("Reference Volume: %.2f m³/ha\n\n", ref_vol_m3_ha))

cat("One-sample t-test: Mean of NS Estimates vs. Reference\n")
print(ns_ttest)

cat("\nOne-sample t-test: Mean of WZ Estimates vs. Reference\n")
print(wz_ttest)

# --- Regression Analysis: Volume vs. Transect Length ---

# 1. Data Preparation for Regression

# According to the grid system, a 10m transect is uniquely identified by 
# the combination of its square number (nr_kwadratu) and line number (nr_linii).
# We combine both NS and WZ datasets and create a full grid of all possible transects
# to ensure empty transects (those with no dead wood) are included in the analysis.

all_transects_raw <- dplyr::bind_rows(
  ns_df %>% dplyr::mutate(dataset = "NS"),
  wz_df %>% dplyr::mutate(dataset = "WZ")
)

# Summarize d2 for transects that have measurements
measured_transects <- all_transects_raw %>%
  dplyr::mutate(d2 = srednica_w_przecieciu^2) %>%
  dplyr::group_by(dataset, nr_kwadratu, nr_linii) %>%
  dplyr::summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop')

# Create a full grid of all possible transects for NS and WZ
ns_grid <- tidyr::expand_grid(dataset = "NS", nr_kwadratu = 1:11, nr_linii = 1:5)
wz_grid <- tidyr::expand_grid(dataset = "WZ", nr_kwadratu = 1:10, nr_linii = 1:6)
full_grid <- dplyr::bind_rows(ns_grid, wz_grid)

# Join measured data with the full grid. Transects without measurements will have sum_d2 = 0.
all_transects <- full_grid %>%
  dplyr::left_join(measured_transects, by = c("dataset", "nr_kwadratu", "nr_linii")) %>%
  dplyr::mutate(sum_d2 = tidyr::replace_na(sum_d2, 0)) %>%
  dplyr::mutate(transect_id = paste(dataset, nr_kwadratu, nr_linii, sep = "_"))


# Function to generate volume estimates for various total transect lengths.
# For lengths greater than a single transect, it combines data from multiple transects.
get_volume_estimates_for_length <- function(total_length_m, transect_pool, transect_unit_length_m) {
  num_transects_to_combine <- total_length_m / transect_unit_length_m
  
  if (num_transects_to_combine > nrow(transect_pool)) {
    return(tibble::tibble(total_length_m = numeric(0), volume_m3_ha = numeric(0)))
  }
  
  if (num_transects_to_combine == 1) {
    result_df <- transect_pool
    result_df$total_length_m <- transect_unit_length_m
    result_df$volume_m3_ha <- (pi^2 * result_df$sum_d2) / (8 * transect_unit_length_m)
    
    return(
      result_df %>%
        dplyr::select("total_length_m", "volume_m3_ha")
    )
  }
  
  n_samples <- 10
  
  if (nrow(transect_pool) >= num_transects_to_combine) {
    purrr::map_df(1:n_samples, ~{
      combined_transects <- transect_pool %>% dplyr::sample_n(num_transects_to_combine, replace = TRUE)
      total_sum_d2 <- sum(combined_transects$sum_d2, na.rm = TRUE)
      
      tibble::tibble(
        total_length_m = total_length_m,
        volume_m3_ha = (pi^2 * total_sum_d2) / (8 * total_length_m)
      )
    })
  } else {
    return(tibble::tibble(total_length_m = numeric(0), volume_m3_ha = numeric(0)))
  }
}

# Generate the regression dataset for transect lengths from 10m up to the maximum possible.
total_transects_available <- nrow(all_transects)
max_length <- total_transects_available * L_PER_TRANSECT

regression_df <- purrr::map_df(
  seq(L_PER_TRANSECT, max_length, by = L_PER_TRANSECT),
  ~get_volume_estimates_for_length(.x, all_transects, L_PER_TRANSECT)
)

# 2. Plotting the Regressions

# --- Linear Model ---
lm_model <- lm(volume_m3_ha ~ total_length_m, data = regression_df)
lm_summary <- summary(lm_model)
p_value_lm <- lm_summary$coefficients[2, 4]
r_squared_lm <- lm_summary$r.squared

# Plot 1: Linear Model with Confidence Interval
lr_plot <- ggplot(regression_df, aes(x = total_length_m, y = volume_m3_ha)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm", aes(color = "Linear Model"), se = TRUE) +
  geom_hline(aes(yintercept = ref_vol_m3_ha, linetype = "Reference"), color = "#004d00") +
  labs(
    title = "Regression of Volume Estimate vs. Transect Length",
    subtitle = "Linear Model with 95% confidence interval",
    x = "Total Transect Length (m)",
    y = "Estimated Volume (m³/ha)"
  ) +
  scale_color_manual(
    name = "Regression Method",
    values = c("Linear Model" = "blue"),
    labels = c(sprintf("Linear Model (R² = %.2f, p = %.3g)", r_squared_lm, p_value_lm))
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Reference" = "dashed"),
    labels = c(sprintf("Reference: %.2f m³/ha", ref_vol_m3_ha))
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave("regression_lm.png", lr_plot, width = 8, height = 6)

# --- LOESS Model ---
# Note: LOESS is non-parametric, so traditional p-values are not applicable.
# We calculate a pseudo R-squared for goodness of fit.
loess_model <- loess(volume_m3_ha ~ total_length_m, data = regression_df)
predictions <- predict(loess_model, newdata = regression_df)
TSS <- sum((regression_df$volume_m3_ha - mean(regression_df$volume_m3_ha))^2)
RSS <- sum((regression_df$volume_m3_ha - predictions)^2, na.rm = TRUE)
pseudo_r_squared_loess <- 1 - (RSS / TSS)

# Plot 2: LOESS with Confidence Interval
loess_plot <- ggplot(regression_df, aes(x = total_length_m, y = volume_m3_ha)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "loess", aes(color = "LOESS"), se = TRUE) +
  geom_hline(aes(yintercept = ref_vol_m3_ha, linetype = "Reference"), color = "#004d00") +
  labs(
    title = "Regression of Volume Estimate vs. Transect Length",
    subtitle = "LOESS smoothing with 95% confidence interval",
    x = "Total Transect Length (m)",
    y = "Estimated Volume (m³/ha)"
  ) +
  scale_color_manual(
    name = "Regression Method",
    values = c("LOESS" = "darkred"),
    labels = c(sprintf("LOESS (Pseudo R² = %.2f)", pseudo_r_squared_loess))
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Reference" = "dashed"),
    labels = c(sprintf("Reference: %.2f m³/ha", ref_vol_m3_ha))
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave("regression_loess.png", loess_plot, width = 8, height = 6)

# --- Quantile Regression ---
# Fit model for the median (tau = 0.5) to get stats
rq_model_median <- quantreg::rq(volume_m3_ha ~ total_length_m, data = regression_df, tau = 0.5)
rq_summary <- summary(rq_model_median, se = "nid") # Use "nid" for p-values
p_value_rq <- rq_summary$coefficients[2, 4]

# Goodness of fit for QR (pseudo R-squared)
null_model <- quantreg::rq(volume_m3_ha ~ 1, tau = 0.5, data = regression_df)
pseudo_r_squared_rq <- 1 - (sum(rq_model_median$rho) / sum(null_model$rho))

# Plot 2: Quantile Regression
quantile_plot <- ggplot(regression_df, aes(x = total_length_m, y = volume_m3_ha)) +
  geom_point(alpha = 0.2) +
  geom_quantile(quantiles = c(0.05, 0.5, 0.95), aes(color = "Quantiles (0.05, 0.5, 0.95)")) +
  geom_hline(aes(yintercept = ref_vol_m3_ha, linetype = "Reference"), color = "#004d00") +
  labs(
    title = "Regression of Volume Estimate vs. Transect Length",
    subtitle = "Uncertainty shown with Quantile Regression",
    x = "Total Transect Length (m)",
    y = "Estimated Volume (m³/ha)"
  ) +
  scale_color_manual(
    name = "Regression Method",
    values = c("Quantiles (0.05, 0.5, 0.95)" = "blue"),
    labels = c(sprintf("Quantiles (0.05, 0.5, 0.95)\nPseudo R² (median) = %.2f, p (median) = %.3g", pseudo_r_squared_rq, p_value_rq))
  ) +
  scale_linetype_manual(
    name = "",
    values = c("Reference" = "dashed"),
    labels = c(sprintf("Reference: %.2f m³/ha", ref_vol_m3_ha))
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical")

ggsave("regression_quantile.png", quantile_plot, width = 8, height = 7)