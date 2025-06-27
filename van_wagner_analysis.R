# Van Wagner Dead Wood Volume Analysis
# Calculate dead wood volume for each group and combined data
# Each group surveyed 500m transect

library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)
library(tidyr)

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

# Read the data
data <- read_csv("data/PELN - transekt liniowy środa.csv")

# Also load NS and WZ data for comparison plots
ns_data <- read_csv("data/PELN - Van Wagner NS (martwe drewno).csv", locale = locale(decimal_mark = ","))
wz_data <- read_csv("data/PELN - Van Wagner WZ (martwe drewno).csv", locale = locale(decimal_mark = ","))

# Clean and prepare the data
# Convert diameter column to numeric, handling comma as decimal separator
data <- data %>%
  mutate(
    Srednica_numeric = case_when(
      `Średnica` == "-" ~ NA_real_,
      `Średnica` == "" ~ NA_real_,
      is.na(`Średnica`) ~ NA_real_,
      TRUE ~ as.numeric(gsub(",", ".", `Średnica`))
    )
  ) %>%
  # Remove rows with missing diameter (class 5 decomposition or unmeasurable)
  filter(!is.na(Srednica_numeric)) %>%
  # Convert diameter from cm to m
  mutate(Srednica_m = Srednica_numeric / 100)

# Van Wagner formula: Volume = π²/8 × Σ(d²)/L
# where d is diameter in meters, L is transect length in meters
van_wagner_constant <- pi^2 / 8
transect_length <- 500  # meters per group

# Calculate volume for each group
group_volumes <- data %>%
  group_by(`Numer grupy (od lewej)`) %>%
  summarise(
    n_logs = n(),
    sum_diameter_squared = sum(Srednica_m^2),
    volume_m3_per_ha = (van_wagner_constant * sum_diameter_squared / transect_length) * 10000,
    .groups = 'drop'
  ) %>%
  arrange(`Numer grupy (od lewej)`)

# Calculate combined volume for all data
total_logs <- nrow(data)
total_groups <- length(unique(data$`Numer grupy (od lewej)`))
total_transect_length <- transect_length * total_groups
total_sum_diameter_squared <- sum(data$Srednica_m^2)
total_volume_m3_per_ha <- (van_wagner_constant * total_sum_diameter_squared / total_transect_length) * 10000

# Display results
cat("VAN WAGNER DEAD WOOD VOLUME ANALYSIS\n")
cat("=====================================\n\n")

cat("Results by Group:\n")
cat("-----------------\n")
for(i in seq_len(nrow(group_volumes))) {
  cat(sprintf("Group %d: %.2f m³/ha (%d logs, transect length: %d m)\n", 
              group_volumes$`Numer grupy (od lewej)`[i], 
              group_volumes$volume_m3_per_ha[i],
              group_volumes$n_logs[i],
              transect_length))
}

cat("\nCombined Results:\n")
cat("------------------\n")
cat(sprintf("Total volume: %.2f m³/ha\n", total_volume_m3_per_ha))
cat(sprintf("Total logs: %d\n", total_logs))
cat(sprintf("Total transect length: %d m (%d groups × %d m)\n", 
            total_transect_length, total_groups, transect_length))
cat(sprintf("Average volume per group: %.2f m³/ha\n", mean(group_volumes$volume_m3_per_ha)))

# Create detailed summary table
summary_table <- group_volumes %>%
  mutate(
    transect_length_m = transect_length,
    volume_formula = paste0("(π²/8) × ", round(sum_diameter_squared, 4), " / ", transect_length, " × 10000")
  ) %>%
  select(
    Group = `Numer grupy (od lewej)`,
    `Number of logs` = n_logs,
    `Sum of d² (m²)` = sum_diameter_squared,
    `Transect length (m)` = transect_length_m,
    `Volume (m³/ha)` = volume_m3_per_ha,
    `Formula` = volume_formula
  )

cat("\nDetailed Summary Table:\n")
cat("=======================\n")
print(summary_table, row.names = FALSE)

# Additional statistics
cat("\nAdditional Statistics:\n")
cat("======================\n")
cat(sprintf("Van Wagner constant (π²/8): %.6f\n", van_wagner_constant))
cat(sprintf("Standard deviation of group volumes: %.2f m³/ha\n", sd(group_volumes$volume_m3_per_ha)))
cat(sprintf("Minimum group volume: %.2f m³/ha (Group %d)\n", 
            min(group_volumes$volume_m3_per_ha),
            group_volumes$`Numer grupy (od lewej)`[which.min(group_volumes$volume_m3_per_ha)]))
cat(sprintf("Maximum group volume: %.2f m³/ha (Group %d)\n", 
            max(group_volumes$volume_m3_per_ha),
            group_volumes$`Numer grupy (od lewej)`[which.max(group_volumes$volume_m3_per_ha)]))

# --- NS and WZ Transect Analysis ---
cat("\nProcessing NS and WZ transect data...\n")

# Constants for NS and WZ analysis
L_PER_TRANSECT_NS <- 50
L_PER_TRANSECT_WZ <- 100

# Debug: Check original column names
cat("Original NS columns:", paste(colnames(ns_data), collapse = ", "), "\n")
cat("Original WZ columns:", paste(colnames(wz_data), collapse = ", "), "\n")

# Process NS data (per square) - use original column name
ns_full_grid <- tibble(nr_kwadratu = 1:11)

ns_estimates <- ns_data %>%
  rename(nr_kwadratu = `Nr kwadratu`) %>%
  # Convert diameter from cm to cm² for Van Wagner formula (diameter is in cm, formula expects cm²)
  mutate(d2 = `Średnica w przecięciu`^2) %>%
  group_by(nr_kwadratu) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  right_join(ns_full_grid, by = "nr_kwadratu") %>%
  mutate(sum_d2 = replace_na(sum_d2, 0)) %>%
  # Van Wagner formula: V(m³/ha) = (π² * Σd²[cm²]) / (8 * L[m])
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT_NS)) %>%
  mutate(method = "NS", transect_id = paste0("Square ", nr_kwadratu))

# Process WZ data (per line) - use original column name
wz_full_grid <- tibble(nr_linii = 1:6)

wz_estimates <- wz_data %>%
  rename(nr_linii = `Nr linii`) %>%
  # Convert diameter from cm to cm² for Van Wagner formula (diameter is in cm, formula expects cm²)
  mutate(d2 = `Średnica w przecięciu`^2) %>%
  group_by(nr_linii) %>%
  summarise(sum_d2 = sum(d2, na.rm = TRUE), .groups = 'drop') %>%
  right_join(wz_full_grid, by = "nr_linii") %>%
  mutate(sum_d2 = replace_na(sum_d2, 0)) %>%
  # Van Wagner formula: V(m³/ha) = (π² * Σd²[cm²]) / (8 * L[m])
  mutate(volume_m3_ha = (pi^2 * sum_d2) / (8 * L_PER_TRANSECT_WZ)) %>%
  mutate(method = "WZ", transect_id = paste0("Line ", nr_linii))

cat(sprintf("NS estimates summary: %.2f ± %.2f m³/ha (n=%d)\n", 
            mean(ns_estimates$volume_m3_ha), sd(ns_estimates$volume_m3_ha), nrow(ns_estimates)))
cat(sprintf("WZ estimates summary: %.2f ± %.2f m³/ha (n=%d)\n", 
            mean(wz_estimates$volume_m3_ha), sd(wz_estimates$volume_m3_ha), nrow(wz_estimates)))

# Save results to CSV
write_csv(summary_table, "van_wagner_results_by_group.csv")

# Create a summary of the combined analysis
combined_summary <- data.frame(
  Analysis = "Combined (All Groups)",
  `Number of logs` = total_logs,
  `Sum of d² (m²)` = total_sum_diameter_squared,
  `Total transect length (m)` = total_transect_length,
  `Volume (m³/ha)` = total_volume_m3_per_ha,
  `Number of groups` = total_groups
)

write_csv(combined_summary, "van_wagner_results_combined.csv")

cat(sprintf("\nResults saved to:\n"))
cat("- van_wagner_results_by_group.csv\n")
cat("- van_wagner_results_combined.csv\n")

# CREATE VISUALIZATIONS
cat("\nCreating visualizations...\n")

# 1. Bar plot of volume by group
p1 <- ggplot(group_volumes, aes(x = factor(`Numer grupy (od lewej)`), y = volume_m3_per_ha)) +
  geom_bar(stat = "identity", fill = BRBG_COLORS$medium_green, color = BRBG_COLORS$dark_green, linewidth = 0.8) +
  geom_text(aes(label = paste0(round(volume_m3_per_ha, 1), "\nm³/ha")), 
            vjust = -0.3, size = 3, fontface = "bold", color = BRBG_COLORS$dark_brown,
            position = position_dodge(width = 0.9)) +
  labs(title = "Dead Wood Volume by Group (Van Wagner Method)",
       subtitle = "Each group surveyed 500m transect",
       x = "Group Number",
       y = "Volume (m³/ha)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = BRBG_COLORS$medium_brown),
        axis.title = element_text(size = 12, color = BRBG_COLORS$dark_brown),
        axis.text = element_text(size = 10, color = BRBG_COLORS$dark_brown),
        panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))  # Add space at top for labels

# Add horizontal line for average
p1 <- p1 + geom_hline(yintercept = mean(group_volumes$volume_m3_per_ha), 
                      linetype = "dashed", color = BRBG_COLORS$medium_brown, linewidth = 1) +
  annotate("text", x = 1, y = mean(group_volumes$volume_m3_per_ha), 
           label = paste0("Avg: ", round(mean(group_volumes$volume_m3_per_ha), 1), " m³/ha"), 
           hjust = 0, vjust = -0.5, color = BRBG_COLORS$medium_brown, size = 3.5, fontface = "bold")

# 1b. Bar plot of volume by NS squares and WZ lines
ns_wz_combined <- bind_rows(
  ns_estimates %>% select(transect_id, volume_m3_ha, method),
  wz_estimates %>% select(transect_id, volume_m3_ha, method)
)

p1b <- ggplot(ns_wz_combined, aes(x = reorder(transect_id, volume_m3_ha), y = volume_m3_ha, fill = method)) +
  geom_bar(stat = "identity", alpha = 0.8, color = "white", linewidth = 0.5) +
  geom_text(aes(label = round(volume_m3_ha, 1)), 
            hjust = -0.1, size = 2.5, fontface = "bold", color = BRBG_COLORS$dark_brown) +
  scale_fill_manual(values = c("NS" = BRBG_COLORS$light_brown, "WZ" = BRBG_COLORS$medium_green)) +
  coord_flip() +
  labs(title = "Dead Wood Volume by NS Squares and WZ Lines",
       subtitle = "Van Wagner method: NS transects 50m, WZ transects 100m",
       x = "Transect",
       y = "Volume (m³/ha)",
       fill = "Method") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = BRBG_COLORS$medium_brown),
        axis.title = element_text(size = 12, color = BRBG_COLORS$dark_brown),
        axis.text = element_text(size = 9, color = BRBG_COLORS$dark_brown),
        legend.position = "bottom",
        legend.title = element_text(color = BRBG_COLORS$dark_brown),
        legend.text = element_text(color = BRBG_COLORS$dark_brown),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

# Add average lines for each method
ns_avg <- mean(ns_estimates$volume_m3_ha)
wz_avg <- mean(wz_estimates$volume_m3_ha)

p1b <- p1b + 
  geom_vline(xintercept = sum(ns_wz_combined$method == "NS") + 0.5, 
             linetype = "dotted", color = "grey60", alpha = 0.7) +
  annotate("text", x = sum(ns_wz_combined$method == "NS")/2, y = max(ns_wz_combined$volume_m3_ha) * 0.9, 
           label = paste0("NS Avg: ", round(ns_avg, 1), " m³/ha"), 
           color = BRBG_COLORS$medium_brown, size = 3, fontface = "bold") +
  annotate("text", x = sum(ns_wz_combined$method == "NS") + sum(ns_wz_combined$method == "WZ")/2, 
           y = max(ns_wz_combined$volume_m3_ha) * 0.9, 
           label = paste0("WZ Avg: ", round(wz_avg, 1), " m³/ha"), 
           color = BRBG_COLORS$medium_brown, size = 3, fontface = "bold")

# 2. Bar plot comparing number of logs vs volume
group_volumes_long <- group_volumes %>%
  select(`Numer grupy (od lewej)`, n_logs, volume_m3_per_ha) %>%
  mutate(
    logs_scaled = n_logs * 2,  # Scale for better visualization
    Group = factor(`Numer grupy (od lewej)`)
  )

p2 <- ggplot(group_volumes_long, aes(x = Group)) +
  geom_bar(aes(y = volume_m3_per_ha, fill = "Volume (m³/ha)"), 
           stat = "identity", alpha = 0.8) +
  geom_bar(aes(y = logs_scaled, fill = "Number of logs (×2)"), 
           stat = "identity", alpha = 0.8) +
  scale_fill_manual(values = c("Volume (m³/ha)" = BRBG_COLORS$medium_green, 
                               "Number of logs (×2)" = BRBG_COLORS$light_brown)) +
  labs(title = "Dead Wood: Volume vs Number of Logs by Group",
       subtitle = "Brown bars show number of logs (scaled ×2 for visualization)",
       x = "Group Number",
       y = "Value",
       fill = "Metric") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = BRBG_COLORS$medium_brown),
        axis.title = element_text(size = 12, color = BRBG_COLORS$dark_brown),
        axis.text = element_text(size = 10, color = BRBG_COLORS$dark_brown),
        legend.position = "bottom",
        legend.title = element_text(color = BRBG_COLORS$dark_brown),
        legend.text = element_text(color = BRBG_COLORS$dark_brown),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3))

# 3. Pie chart of volume contribution by group
# Create BrBG palette for groups
group_colors <- c(BRBG_COLORS$dark_brown, BRBG_COLORS$medium_brown, 
                  BRBG_COLORS$light_brown, BRBG_COLORS$medium_green,
                  BRBG_COLORS$dark_green)

p3 <- ggplot(group_volumes, aes(x = "", y = volume_m3_per_ha, 
                                fill = factor(`Numer grupy (od lewej)`))) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 2) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Proportion of Dead Wood Volume by Group",
       fill = "Group") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        legend.title = element_text(color = BRBG_COLORS$dark_brown),
        legend.text = element_text(color = BRBG_COLORS$dark_brown),
        plot.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(20, 20, 20, 20, "pt")) +
  geom_text(aes(label = paste0("Group ", `Numer grupy (od lewej)`, "\n", 
                               round(volume_m3_per_ha, 1), " m³/ha")), 
            position = position_stack(vjust = 0.5), size = 3, 
            color = "white", fontface = "bold")

# 4. Summary statistics plot
stats_df <- data.frame(
  Metric = c("Total Volume", "Average per Group", "Standard Deviation", "Range"),
  Value = c(total_volume_m3_per_ha, 
            mean(group_volumes$volume_m3_per_ha),
            sd(group_volumes$volume_m3_per_ha),
            max(group_volumes$volume_m3_per_ha) - min(group_volumes$volume_m3_per_ha)),
  Unit = c("m³/ha", "m³/ha", "m³/ha", "m³/ha")
)

p4 <- ggplot(stats_df, aes(x = reorder(Metric, Value), y = Value)) +
  geom_bar(stat = "identity", fill = BRBG_COLORS$dark_green, alpha = 0.8, 
           color = BRBG_COLORS$medium_green, linewidth = 0.8) +
  geom_text(aes(label = paste0(round(Value, 1), " ", Unit)), 
            hjust = -0.1, size = 3.5, fontface = "bold", color = BRBG_COLORS$dark_brown) +
  coord_flip() +
  labs(title = "Summary Statistics - Dead Wood Volume",
       x = "Metric",
       y = "Value (m³/ha)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        axis.title = element_text(size = 12, color = BRBG_COLORS$dark_brown),
        axis.text = element_text(size = 10, color = BRBG_COLORS$dark_brown),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
        panel.grid.minor = element_line(color = "grey95", linewidth = 0.3))

# 5. Detailed comparison table plot
table_data <- group_volumes %>%
  mutate(
    Group = paste("Group", `Numer grupy (od lewej)`),
    Volume = paste0(round(volume_m3_per_ha, 2), " m³/ha"),
    Logs = paste(n_logs, "logs"),
    `Sum d²` = paste0(round(sum_diameter_squared, 3), " m²")
  ) %>%
  select(Group, Volume, Logs, `Sum d²`)

# Create table plot
library(gridExtra)
p5 <- tableGrob(table_data, rows = NULL, 
                theme = ttheme_default(base_size = 10,
                                     core = list(bg_params = list(fill = BRBG_COLORS$light_gray),
                                                fg_params = list(col = BRBG_COLORS$dark_brown)),
                                     colhead = list(bg_params = list(fill = BRBG_COLORS$medium_green),
                                                   fg_params = list(col = "white", fontface = "bold"))))

# Save individual plots
ggsave("van_wagner_volume_by_group.png", p1, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("van_wagner_ns_wz_transects.png", p1b, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("van_wagner_volume_vs_logs.png", p2, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("van_wagner_volume_pie.png", p3, width = 8, height = 8, dpi = 300, bg = "white")
ggsave("van_wagner_summary_stats.png", p4, width = 10, height = 6, dpi = 300, bg = "white")

# Create combined plot
combined_plot <- grid.arrange(
  p1, p1b,
  p3, p4,
  ncol = 2,
  top = textGrob("Van Wagner Dead Wood Volume Analysis - Complete Results", 
                 gp = gpar(fontsize = 16, fontface = "bold", col = BRBG_COLORS$dark_brown))
)

ggsave("van_wagner_combined_analysis.png", combined_plot, width = 16, height = 12, dpi = 300, bg = "white")

# Save table as separate plot
png("van_wagner_results_table.png", width = 800, height = 300, res = 150, bg = "white")
grid.draw(p5)
dev.off()

cat("\nPlots saved:\n")
cat("- van_wagner_volume_by_group.png\n")
cat("- van_wagner_ns_wz_transects.png\n")
cat("- van_wagner_volume_vs_logs.png\n")
cat("- van_wagner_volume_pie.png\n")
cat("- van_wagner_summary_stats.png\n")
cat("- van_wagner_combined_analysis.png\n")
cat("- van_wagner_results_table.png\n")

# Display final summary
cat("\n", rep("=", 50), "\n", sep = "")
cat("FINAL SUMMARY\n")
cat(rep("=", 50), "\n", sep = "")
cat(sprintf("Total dead wood volume: %.2f m³/ha\n", total_volume_m3_per_ha))
cat(sprintf("Volume range: %.2f - %.2f m³/ha\n", 
            min(group_volumes$volume_m3_per_ha), 
            max(group_volumes$volume_m3_per_ha)))
cat(sprintf("Most productive group: Group %d (%.2f m³/ha)\n", 
            group_volumes$`Numer grupy (od lewej)`[which.max(group_volumes$volume_m3_per_ha)],
            max(group_volumes$volume_m3_per_ha)))
cat(sprintf("Least productive group: Group %d (%.2f m³/ha)\n", 
            group_volumes$`Numer grupy (od lewej)`[which.min(group_volumes$volume_m3_per_ha)],
            min(group_volumes$volume_m3_per_ha)))
cat(rep("=", 50), "\n", sep = "")
