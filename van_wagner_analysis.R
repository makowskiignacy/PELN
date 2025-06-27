# Van Wagner Dead Wood Volume Analysis
# Calculate dead wood volume for each group and combined data
# Each group surveyed 500m transect

library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(grid)
library(RColorBrewer)

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
  geom_bar(stat = "identity", fill = BRBG_COLORS$medium_green, alpha = 0.8, color = BRBG_COLORS$dark_green, linewidth = 0.8) +
  geom_text(aes(label = paste0(round(volume_m3_per_ha, 1), "\nm³/ha")), 
            vjust = -0.5, size = 3.5, fontface = "bold", color = BRBG_COLORS$dark_brown) +
  labs(title = "Dead Wood Volume by Group (Van Wagner Method)",
       subtitle = "Each group surveyed 500m transect",
       x = "Group Number",
       y = "Volume (m³/ha)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", color = BRBG_COLORS$dark_brown),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = BRBG_COLORS$medium_brown),
        axis.title = element_text(size = 12, color = BRBG_COLORS$dark_brown),
        axis.text = element_text(size = 10, color = BRBG_COLORS$dark_brown),
        panel.background = element_rect(fill = BRBG_COLORS$light_gray, color = NA))

# Add horizontal line for average
p1 <- p1 + geom_hline(yintercept = mean(group_volumes$volume_m3_per_ha), 
                      linetype = "dashed", color = BRBG_COLORS$medium_brown, alpha = 0.8, linewidth = 1) +
  annotate("text", x = Inf, y = mean(group_volumes$volume_m3_per_ha), 
           label = paste0("Average: ", round(mean(group_volumes$volume_m3_per_ha), 1), " m³/ha"), 
           hjust = 1.1, vjust = -0.5, color = BRBG_COLORS$medium_brown, size = 3, fontface = "bold")

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
        panel.background = element_rect(fill = BRBG_COLORS$light_gray, color = NA))

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
        legend.text = element_text(color = BRBG_COLORS$dark_brown)) +
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
        panel.background = element_rect(fill = BRBG_COLORS$light_gray, color = NA))

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
                                     core = list(bg_params = list(fill = BRBG_COLORS$light_gray, alpha = 0.8),
                                                fg_params = list(col = BRBG_COLORS$dark_brown)),
                                     colhead = list(bg_params = list(fill = BRBG_COLORS$medium_green, alpha = 0.9),
                                                   fg_params = list(col = "white", fontface = "bold"))))

# Save individual plots
ggsave("van_wagner_volume_by_group.png", p1, width = 10, height = 6, dpi = 300)
ggsave("van_wagner_volume_vs_logs.png", p2, width = 10, height = 6, dpi = 300)
ggsave("van_wagner_volume_pie.png", p3, width = 8, height = 8, dpi = 300)
ggsave("van_wagner_summary_stats.png", p4, width = 10, height = 6, dpi = 300)

# Create combined plot
combined_plot <- grid.arrange(
  p1, p2,
  p3, p4,
  ncol = 2,
  top = textGrob("Van Wagner Dead Wood Volume Analysis - Complete Results", 
                 gp = gpar(fontsize = 16, fontface = "bold", col = BRBG_COLORS$dark_brown))
)

ggsave("van_wagner_combined_analysis.png", combined_plot, width = 16, height = 12, dpi = 300)

# Save table as separate plot
png("van_wagner_results_table.png", width = 800, height = 300, res = 150)
grid.draw(p5)
dev.off()

cat("\nPlots saved:\n")
cat("- van_wagner_volume_by_group.png\n")
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
