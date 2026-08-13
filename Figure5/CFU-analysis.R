# ============================================================
# Figure 5 — CFU analysis and plots for T0 and T1
# ============================================================

# ---- Libraries ----
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(multcompView)
library(car)
library(paletteer)

# ---- Settings ----
input_file <- "data/time1-repeat-cfu.csv"
output_dir <- "figures_test"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

samples <- c(
  "W16.3_WT",
  "W16.3_GacA?",
  "W16.3_GacS?",
  "WT:GacA_1:2",
  "WT:GacA_2:1",
  "W15.13_WT",
  "W23.13_WT"
)

sample_order <- c(
  "W16.3_WT",
  "W23.13_WT",
  "W15.13_WT",
  "W16.3_GacA?",
  "W16.3_GacS?",
  "WT:GacA_2:1",
  "WT:GacA_1:2"
)

sample_labels <- c(
  "W16.3 WT",
  "W23.13 WT",
  "W15.13 WT",
  "W16.3 GacA∆",
  "W16.3 GacS∆",
  "W16.3 WT:GacA∆ 2:1",
  "W16.3 WT:GacA∆ 1:2"
)

# ---- Read data ----
df <- read_csv(input_file, show_col_types = FALSE)


# ============================================================
# Function to analyse one time point
# ============================================================

analyse_timepoint <- function(data, timepoint) {
  
  message("Analysing ", timepoint, "...")
  
  # ----------------------------------------------------------
  # Prepare data
  # ----------------------------------------------------------
  
  df_mean <- data %>%
    filter(
      Sample %in% samples,
      time == timepoint
    ) %>%
    mutate(
      log_cfu_ml = as.numeric(log_cfu_ml),
      Replicate = factor(Replicate),
      Sample = factor(Sample, levels = sample_order)
    ) %>%
    drop_na(log_cfu_ml) %>%
    group_by(
      Sample,
      Tree,
      Inoculation,
      Replicate
    ) %>%
    summarise(
      mean_log_cfu_ml = mean(log_cfu_ml, na.rm = TRUE),
      sd = sd(log_cfu_ml, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # Descriptive statistics
  # ----------------------------------------------------------
  
  descriptive_stats <- df_mean %>%
    group_by(Sample) %>%
    summarise(
      n = n(),
      mean = mean(mean_log_cfu_ml, na.rm = TRUE),
      sd = sd(mean_log_cfu_ml, na.rm = TRUE),
      median = median(mean_log_cfu_ml, na.rm = TRUE),
      min = min(mean_log_cfu_ml, na.rm = TRUE),
      max = max(mean_log_cfu_ml, na.rm = TRUE),
      se = sd / sqrt(n),
      .groups = "drop"
    )
  
  
  # ----------------------------------------------------------
  # ANOVA
  # ----------------------------------------------------------
  
  anova_model <- aov(
    mean_log_cfu_ml ~ Sample,
    data = df_mean
  )
  
  anova_summary <- summary(anova_model)[[1]]
  
  anova_p <- anova_summary[["Pr(>F)"]][1]
  
  
  # ----------------------------------------------------------
  # Assumption checks
  # ----------------------------------------------------------
  
  residuals_anova <- residuals(anova_model)
  
  shapiro_result <- shapiro.test(residuals_anova)
  
  levene_result <- car::leveneTest(
    mean_log_cfu_ml ~ Sample,
    data = df_mean
  )
  
  assumption_stats <- tibble(
    time = timepoint,
    shapiro_W = unname(shapiro_result$statistic),
    shapiro_p = shapiro_result$p.value,
    levene_F = levene_result$`F value`[1],
    levene_p = levene_result$`Pr(>F)`[1]
  )
  
  
  # ----------------------------------------------------------
  # Tukey HSD
  # ----------------------------------------------------------
  
  tukey <- TukeyHSD(anova_model)
  
  tukey_df <- as.data.frame(tukey$Sample) %>%
    tibble::rownames_to_column("Comparison") %>%
    rename(
      diff = diff,
      lwr = lwr,
      upr = upr,
      p_adj = `p adj`
    ) %>%
    mutate(
      time = timepoint,
      significant = p_adj <= 0.05
    ) %>%
    select(
      time,
      Comparison,
      diff,
      lwr,
      upr,
      p_adj,
      significant
    )
  
  
  # ----------------------------------------------------------
  # Compact letter display
  # ----------------------------------------------------------
  
  cld <- multcompLetters4(
    anova_model,
    tukey
  )
  
  cld_letters <- tibble(
    Sample = names(cld$Sample$Letters),
    cld = unname(cld$Sample$Letters)
  )
  
  label_positions <- df_mean %>%
    group_by(Sample) %>%
    summarise(
      y = max(mean_log_cfu_ml, na.rm = TRUE) + 0.3,
      .groups = "drop"
    )
  
  cld_letters <- cld_letters %>%
    left_join(
      label_positions,
      by = "Sample"
    ) %>%
    mutate(
      Sample = factor(
        Sample,
        levels = sample_order
      )
    )
  
  
  # ----------------------------------------------------------
  # ANOVA statistics table
  # ----------------------------------------------------------
  
  anova_stats <- tibble(
    time = timepoint,
    term = "Sample",
    df = anova_summary[["Df"]][1],
    sum_sq = anova_summary[["Sum Sq"]][1],
    mean_sq = anova_summary[["Mean Sq"]][1],
    F_value = anova_summary[["F value"]][1],
    p_value = anova_p
  )
  
  
  # ----------------------------------------------------------
  # Combined statistics table
  # ----------------------------------------------------------
  
  stats_table <- descriptive_stats %>%
    left_join(
      cld_letters %>%
        select(Sample, cld),
      by = "Sample"
    ) %>%
    mutate(
      time = timepoint
    ) %>%
    select(
      time,
      Sample,
      n,
      mean,
      sd,
      median,
      min,
      max,
      se,
      cld
    )
  
  
  # ----------------------------------------------------------
  # Save statistics
  # ----------------------------------------------------------
  
  write_csv(
    stats_table,
    file.path(
      output_dir,
      paste0(
        "Figure5-stats-",
        tolower(timepoint),
        ".csv"
      )
    )
  )
  
  write_csv(
    tukey_df,
    file.path(
      output_dir,
      paste0(
        "Figure5-tukey-",
        tolower(timepoint),
        ".csv"
      )
    )
  )
  
  write_csv(
    anova_stats,
    file.path(
      output_dir,
      paste0(
        "Figure5-anova-",
        tolower(timepoint),
        ".csv"
      )
    )
  )
  
  write_csv(
    assumption_stats,
    file.path(
      output_dir,
      paste0(
        "Figure5-assumptions-",
        tolower(timepoint),
        ".csv"
      )
    )
  )
  
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  y_min <- min(
    4,
    min(
      df_mean$mean_log_cfu_ml,
      na.rm = TRUE
    )
  )
  
  y_max <- max(
    df_mean$mean_log_cfu_ml,
    na.rm = TRUE
  ) + 0.7
  
  
  p <- ggplot(
    df_mean,
    aes(
      x = Sample,
      y = mean_log_cfu_ml
    )
  ) +
    
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65
    ) +
    
    geom_jitter(
      aes(colour = Replicate),
      width = 0.10,
      size = 1.8,
      alpha = 0.9
    ) +
    
    geom_text(
      data = cld_letters,
      aes(
        x = Sample,
        y = y,
        label = cld
      ),
      inherit.aes = FALSE,
      size = 3.5
    ) +
    
    theme_classic() +
    
    scale_colour_paletteer_d(
      "tvthemes::Alexandrite"
    ) +
    
    scale_x_discrete(
      labels = sample_labels
    ) +
    
    coord_cartesian(
      ylim = c(y_min, y_max)
    ) +
    
    labs(
      title = paste0(
        "Log CFU/ml — ",
        timepoint
      ),
      x = "Sample",
      y = "Log (CFU ml⁻¹ + 1)",
      colour = "Replicate"
    ) +
    
    theme(
      axis.text.x = element_text(
        size = 10,
        vjust = 0.5,
        hjust = 1,
        angle = 90
      ),
      axis.title = element_text(
        size = 10
      ),
      axis.text.y = element_text(
        size = 10
      ),
      plot.title = element_text(
        size = 11
      ),
      legend.text = element_text(
        size = 10
      ),
      legend.title = element_text(
        size = 10
      )
    )
  
  
  # ----------------------------------------------------------
  # Save plot
  # ----------------------------------------------------------
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0(
        "Figure5-",
        tolower(timepoint),
        ".png"
      )
    ),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  
  # ----------------------------------------------------------
  # Print results
  # ----------------------------------------------------------
  
  print(p)
  
  message(
    "\n",
    timepoint,
    " ANOVA p = ",
    signif(anova_p, 4),
    "\nShapiro-Wilk p = ",
    signif(shapiro_result$p.value, 4),
    "\nLevene's test p = ",
    signif(
      levene_result$`Pr(>F)`[1],
      4
    ),
    "\n"
  )
  
  invisible(
    list(
      data = df_mean,
      descriptive = stats_table,
      anova = anova_stats,
      tukey = tukey_df,
      assumptions = assumption_stats,
      cld = cld_letters,
      plot = p
    )
  )
}


# ============================================================
# Run analysis for BOTH time points
# ============================================================

results_t0 <- analyse_timepoint(
  df,
  "T0"
)

results_t1 <- analyse_timepoint(
  df,
  "T1"
)

message("Done.")
message(
  "Statistics and plots saved in: ",
  normalizePath(output_dir)
)

