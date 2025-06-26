library(tidyverse)
library(janitor)

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