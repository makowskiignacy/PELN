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
L_PER_TRANSECT_NS <- 50
L_PER_TRANSECT_WZ <- 100

# --- BrBG Color Palette Definition ---
# Define consistent BrBG colors for all plots
BRBG_COLORS <- list(
  dark_brown = "#543005",
  medium_brown = "#8c510a", 
  light_brown = "#bf812d",
  dark_green = "#003c30",
  medium_green = "#35978f",
  light_gray = "#f5f5f5"
)

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
# Create full grid of all possible squares (1-11) to include empty transects
ns_full_grid <- tibble(nr_kwadratu = 1:11)

ns_estimates <- ns_df %>%
  mutate(d2 = srednica_w_przecieciu^2) %>%
  group_by(nr_kwadratu) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  right_join(ns_full_grid, by = "nr_kwadratu") %>%
  mutate(sum_d2 = replace_na(sum_d2, 0)) %>%
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT_NS))

# --- WZ Volume Estimation (per Line) ---
# Create full grid of all possible lines (1-6) to include empty transects
wz_full_grid <- tibble(nr_linii = 1:6)

wz_estimates <- wz_df %>%
  mutate(d2 = srednica_w_przecieciu^2) %>%
  group_by(nr_linii) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  right_join(wz_full_grid, by = "nr_linii") %>%
  mutate(sum_d2 = replace_na(sum_d2, 0)) %>%
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT_WZ))

# --- Data Preparation for Plotting ---
plot_df <- bind_rows(
  ns_estimates %>% select(volume_m3_ha) %>% mutate(method = "NS per Square"),
  wz_estimates %>% select(volume_m3_ha) %>% mutate(method = "WZ per Line")
)

# --- DEBUG: Print all transect lengths ---
cat("--- DEBUG: All transect volume estimates ---\n")
cat("NS estimates (per square):\n")
print(ns_estimates %>% select(nr_kwadratu, volume_m3_ha) %>% arrange(nr_kwadratu))
cat("\nWZ estimates (per line):\n") 
print(wz_estimates %>% select(nr_linii, volume_m3_ha) %>% arrange(nr_linii))
cat("\n")

# --- Plotting ---
# Define BrBG color palette using consistent colors
elegant_palette <- c("NS per Square" = BRBG_COLORS$dark_brown, "WZ per Line" = BRBG_COLORS$medium_brown)

violin_plot <- ggplot(plot_df, aes(x = method, y = volume_m3_ha)) +
  stat_ydensity(geom = "violin", aes(fill = method), trim = FALSE, alpha = 0.6, bounds = c(0, Inf)) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.5) +
  # Add a point marker for the reference value on each violin plot
  geom_point(data = . %>% distinct(method), aes(y = ref_vol_m3_ha, shape = "Reference"), size = 5, color = BRBG_COLORS$dark_green) +
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
  geom_point(alpha = 0.2, color = BRBG_COLORS$light_brown) +
  geom_quantile(quantiles = c(0.05, 0.5, 0.95), aes(color = "Quantiles (0.05, 0.5, 0.95)")) +
  geom_hline(aes(yintercept = ref_vol_m3_ha, linetype = "Reference"), color = BRBG_COLORS$dark_green, linewidth = 1.3) +
  labs(
    title = "Regression of Volume Estimate vs. Transect Length",
    subtitle = "Uncertainty shown with Quantile Regression",
    x = "Total Transect Length (m)",
    y = "Estimated Volume (m³/ha)"
  ) +
  scale_color_manual(
    name = "Regression Method",
    values = c("Quantiles (0.05, 0.5, 0.95)" = BRBG_COLORS$dark_brown),
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


# --- Cumulative Volume Analysis ---

# This analysis calculates the dead wood volume estimate cumulatively,
# starting from the first transect and progressively adding more.
# It helps visualize how the estimate stabilizes as more data is included.

# Order transects for a consistent cumulative calculation
ns_transects_ordered <- all_transects %>%
  dplyr::filter(dataset == "NS") %>%
  dplyr::arrange(nr_kwadratu, nr_linii)

wz_transects_ordered <- all_transects %>%
  dplyr::filter(dataset == "WZ") %>%
  dplyr::arrange(nr_kwadratu, nr_linii)

# Calculate cumulative volume for NS
ns_cumulative_df <- ns_transects_ordered %>%
  dplyr::mutate(
    cumulative_sum_d2 = cumsum(sum_d2),
    transect_count = dplyr::row_number(),
    total_length_m = transect_count * L_PER_TRANSECT,
    volume_m3_ha = (pi^2 * cumulative_sum_d2) / (8 * total_length_m)
  ) %>%
  dplyr::mutate(method = "NS")

# Calculate cumulative volume for WZ
wz_cumulative_df <- wz_transects_ordered %>%
  dplyr::mutate(
    cumulative_sum_d2 = cumsum(sum_d2),
    transect_count = dplyr::row_number(),
    total_length_m = transect_count * L_PER_TRANSECT,
    volume_m3_ha = (pi^2 * cumulative_sum_d2) / (8 * total_length_m)
  ) %>%
  dplyr::mutate(method = "WZ")

# --- Combined Cumulative Volume and Uncertainty ---

# Create alternating pattern between NS and WZ transects
ns_transects_for_alternating <- all_transects %>%
  dplyr::filter(dataset == "NS") %>%
  dplyr::arrange(nr_kwadratu, nr_linii)

wz_transects_for_alternating <- all_transects %>%
  dplyr::filter(dataset == "WZ") %>%
  dplyr::arrange(nr_kwadratu, nr_linii)

# Create alternating combined dataset
create_alternating_combined_transects <- function() {
  max_length <- max(nrow(ns_transects_for_alternating), nrow(wz_transects_for_alternating))
  combined_alternating <- tibble()
  
  for (i in 1:max_length) {
    # Add NS transect if available
    if (i <= nrow(ns_transects_for_alternating)) {
      combined_alternating <- bind_rows(combined_alternating, ns_transects_for_alternating[i, ])
    }
    # Add WZ transect if available
    if (i <= nrow(wz_transects_for_alternating)) {
      combined_alternating <- bind_rows(combined_alternating, wz_transects_for_alternating[i, ])
    }
  }
  return(combined_alternating)
}

alternating_transects <- create_alternating_combined_transects()

# Calculate cumulative volume for alternating combined dataset
combined_cumulative_df <- alternating_transects %>%
  dplyr::mutate(
    cumulative_sum_d2 = cumsum(sum_d2),
    transect_count = dplyr::row_number(),
    total_length_m = transect_count * L_PER_TRANSECT,
    volume_m3_ha = (pi^2 * cumulative_sum_d2) / (8 * total_length_m)
  ) %>%
  dplyr::mutate(method = "Combined")

# --- Bootstrap for Uncertainty of Combined Estimate ---
# We simulate different accumulation paths by creating different alternating patterns
n_bootstrap <- 500 # Number of bootstrap samples.

# Function to get one bootstrap sample with different alternating pattern
get_bootstrap_alternating_volume <- function(seed) {
  set.seed(seed)
  
  # Shuffle both datasets independently
  ns_shuffled <- ns_transects_for_alternating[sample(seq_len(nrow(ns_transects_for_alternating))), ]
  wz_shuffled <- wz_transects_for_alternating[sample(seq_len(nrow(wz_transects_for_alternating))), ]
  
  # Create alternating pattern with shuffled data
  max_length <- max(nrow(ns_shuffled), nrow(wz_shuffled))
  bootstrap_combined <- tibble()
  
  for (i in 1:max_length) {
    if (i <= nrow(ns_shuffled)) {
      bootstrap_combined <- bind_rows(bootstrap_combined, ns_shuffled[i, ])
    }
    if (i <= nrow(wz_shuffled)) {
      bootstrap_combined <- bind_rows(bootstrap_combined, wz_shuffled[i, ])
    }
  }
  
  # Calculate cumulative volume
  bootstrap_combined %>%
    dplyr::mutate(
      cumulative_sum_d2 = cumsum(sum_d2),
      transect_count = dplyr::row_number(),
      total_length_m = transect_count * L_PER_TRANSECT,
      volume_m3_ha = (pi^2 * cumulative_sum_d2) / (8 * total_length_m)
    ) %>%
    dplyr::select(total_length_m, volume_m3_ha) %>%
    dplyr::mutate(bootstrap_run = seed)
}

# Generate bootstrap replicates
bootstrap_runs <- purrr::map_df(1:n_bootstrap, ~get_bootstrap_alternating_volume(.x))

# Calculate confidence intervals and relative uncertainty from bootstrap runs
cumulative_ci <- bootstrap_runs %>%
  dplyr::group_by(total_length_m) %>%
  dplyr::summarise(
    mean_volume = mean(volume_m3_ha, na.rm = TRUE),
    lower_ci = quantile(volume_m3_ha, 0.025, na.rm = TRUE),
    upper_ci = quantile(volume_m3_ha, 0.975, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  dplyr::mutate(
    ci_width = upper_ci - lower_ci,
    relative_uncertainty = (ci_width / (2 * mean_volume)) * 100  # Half-width as percentage of mean
  )

# Find where uncertainty drops below 10%
uncertainty_threshold <- cumulative_ci %>%
  dplyr::filter(relative_uncertainty <= 10) %>%
  dplyr::slice(1)

# Join CIs with the main combined cumulative dataframe
combined_cumulative_df_with_ci <- combined_cumulative_df %>%
  dplyr::left_join(cumulative_ci, by = "total_length_m")


# Combine for plotting
cumulative_plot_df <- dplyr::bind_rows(ns_cumulative_df, wz_cumulative_df)

# Create segments for the uncertainty heatmap
create_uncertainty_segments <- function(data) {
  segments <- data.frame()
  for(i in 1:(nrow(data)-1)) {
    segments <- rbind(segments, data.frame(
      x_start = data$total_length_m[i],
      x_end = data$total_length_m[i+1],
      y_lower = data$lower_ci[i],
      y_upper = data$upper_ci[i],
      uncertainty = data$relative_uncertainty[i]
    ))
  }
  return(segments)
}

uncertainty_segments <- create_uncertainty_segments(combined_cumulative_df_with_ci)

# Create the plot with uncertainty heatmap using segments
cumulative_volume_plot <- ggplot() +
  # Plot uncertainty as colored segments
  geom_rect(data = uncertainty_segments, 
            aes(xmin = x_start, xmax = x_end, ymin = y_lower, ymax = y_upper, fill = uncertainty),
            alpha = 0.4) +
  
  # Plot NS and WZ lines
  geom_line(data = cumulative_plot_df, aes(x = total_length_m, y = volume_m3_ha, color = method), linewidth = 1.1) +
  geom_point(data = cumulative_plot_df, aes(x = total_length_m, y = volume_m3_ha, color = method), alpha = 0.6) +
  
  # Plot Combined line
  geom_line(data = combined_cumulative_df_with_ci, aes(x = total_length_m, y = volume_m3_ha, color = "Combined"), linewidth = 1.2) +
  
  # Mark the point where uncertainty drops below 10%
  {if(nrow(uncertainty_threshold) > 0) {
    geom_vline(xintercept = uncertainty_threshold$total_length_m, linetype = "dotted", color = BRBG_COLORS$medium_brown, linewidth = 1)
  }} +
  {if(nrow(uncertainty_threshold) > 0) {
    geom_point(data = uncertainty_threshold, aes(x = total_length_m, y = mean_volume), 
               color = BRBG_COLORS$medium_brown, size = 4, shape = 17)
  }} +
  
  # Plot reference line
  geom_hline(aes(yintercept = ref_vol_m3_ha, linetype = "Reference"), color = BRBG_COLORS$dark_green, linewidth = 1.2) +
  
  labs(
    title = "Cumulative Dead Wood Volume Estimate vs. Transect Length",
    subtitle = paste0("Alternating NS-WZ sampling pattern with uncertainty heatmap",
                     if(nrow(uncertainty_threshold) > 0) 
                       paste0("\nUncertainty ≤10% achieved at ", uncertainty_threshold$total_length_m, "m") 
                     else "\nUncertainty never drops below 10%"),
    x = "Total Transect Length (m)",
    y = "Estimated Volume (m³/ha)",
    color = "Method",
    fill = "Relative Uncertainty (%)"
  ) +
  scale_color_manual(name = "Method", values = c("NS" = BRBG_COLORS$dark_brown, "WZ" = BRBG_COLORS$medium_brown, "Combined" = BRBG_COLORS$dark_green)) +
  scale_fill_gradient2(name = "Uncertainty (%)", 
                       low = BRBG_COLORS$medium_green, mid = BRBG_COLORS$light_gray, high = BRBG_COLORS$light_brown, 
                       midpoint = 15, limits = c(0, NA),
                       guide = guide_colorbar(title.position = "top")) +
  scale_linetype_manual(
    name = "",
    values = c("Reference" = "dashed"),
    labels = c(sprintf("Reference: %.2f m³/ha", ref_vol_m3_ha))
  ) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "vertical")

# Print summary statistics
if(nrow(uncertainty_threshold) > 0) {
  cat(sprintf("\n--- Uncertainty Analysis ---\n"))
  cat(sprintf("Uncertainty drops below 10%% at transect length: %.0f m\n", uncertainty_threshold$total_length_m))
  cat(sprintf("Volume estimate at this point: %.2f m³/ha (±%.1f%%)\n", 
              uncertainty_threshold$mean_volume, uncertainty_threshold$relative_uncertainty))
} else {
  cat("\n--- Uncertainty Analysis ---\n")
  cat("Uncertainty never drops below 10% within the available transect length.\n")
  min_uncertainty <- min(cumulative_ci$relative_uncertainty, na.rm = TRUE)
  cat(sprintf("Minimum uncertainty achieved: %.1f%%\n", min_uncertainty))
}


# Save the plot
ggsave("cumulative_volume_plot.png", cumulative_volume_plot, width = 10, height = 7)