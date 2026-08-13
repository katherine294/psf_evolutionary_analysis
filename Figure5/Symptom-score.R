# ============================================================
# Figure 5 — Symptom analysis at T1
# ============================================================

# ---- Libraries ----
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)
library(viridis)


# ---- Settings ----
input_file <- "data/cfu_and_symptoms.csv"
output_dir <- "symptom-figures"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ---- Sample order ----
sample_order <- c(
  "PBS",
  "W16.3_WT",
  "W23.13_WT",
  "W15.13_WT",
  "W16.3_GacA?",
  "W16.3_GacS?",
  "WT:GacA_1:2",
  "WT:GacA_2:1"
)


# ---- Plot labels ----
sample_labels <- c(
  "PBS",
  "W16.3 WT",
  "W23.13 WT",
  "W15.13 WT",
  "W16.3 GacA∆",
  "W16.3 GacS∆",
  "W16.3 WT:GacA∆ 1:2",
  "W16.3 WT:GacA∆ 2:1"
)


# ============================================================
# Read and prepare data
# ============================================================

df <- read_csv(
  input_file,
  show_col_types = FALSE
) %>%
  filter(
      Sample %in% sample_order
  ) %>%
  mutate(
    Sample = factor(
      Sample,
      levels = sample_order
    ),
    Symptom_0_3 = as.numeric(Symptom_0_3)
  ) %>%
  drop_na(
    Sample,
    Symptom_0_3
  )


# ============================================================
# Kruskal-Wallis test
# ============================================================

kruskal_result <- kruskal.test(
  Symptom_0_3 ~ Sample,
  data = df
)

print(kruskal_result)


# ============================================================
# Dunn's post-hoc test
# Benjamini-Hochberg correction
# ============================================================

dunn_results <- dunnTest(
  Symptom_0_3 ~ Sample,
  data = df,
  method = "bh"
)


# Extract results
dunn_df <- dunn_results$res %>%
  mutate(
    significant = P.adj <= 0.05
  ) %>%
  select(
    Comparison,
    Z,
    P.unadj,
    P.adj,
    significant
  )


# Significant comparisons only
significant_dunn <- dunn_df %>%
  filter(
    significant
  )


# ============================================================
# Save Dunn test results
# ============================================================

write_csv(
  dunn_df,
  file.path(
    output_dir,
    "Figure5-symptom-dunn-t1.csv"
  )
)

write_csv(
  significant_dunn,
  file.path(
    output_dir,
    "Figure5-symptom-dunn-significant-t1.csv"
  )
)


# ============================================================
# Save Kruskal-Wallis result
# ============================================================

kruskal_stats <- tibble(
  test = "Kruskal-Wallis",
  statistic = unname(
    kruskal_result$statistic
  ),
  df = unname(
    kruskal_result$parameter
  ),
  p_value = kruskal_result$p.value
)

write_csv(
  kruskal_stats,
  file.path(
    output_dir,
    "Figure5-symptom-kruskal-t1.csv"
  )
)


# ============================================================
# Symptom score counts
# ============================================================

score_counts <- df %>%
  count(
    Sample,
    Symptom_0_3
  ) %>%
  complete(
    Sample,
    Symptom_0_3 = 0:3,
    fill = list(n = 0)
  )


# ============================================================
# Descriptive statistics
# ============================================================

symptom_stats <- df %>%
  group_by(Sample) %>%
  summarise(
    n = n(),
    mean = mean(
      Symptom_0_3,
      na.rm = TRUE
    ),
    sd = sd(
      Symptom_0_3,
      na.rm = TRUE
    ),
    median = median(
      Symptom_0_3,
      na.rm = TRUE
    ),
    min = min(
      Symptom_0_3,
      na.rm = TRUE
    ),
    max = max(
      Symptom_0_3,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

write_csv(
  symptom_stats,
  file.path(
    output_dir,
    "Figure5-symptom-stats-t1.csv"
  )
)


# ============================================================
# Symptom graph
# ============================================================

p_symptoms <- ggplot(
  score_counts,
  aes(
    x = Sample,
    y = n,
    fill = factor(Symptom_0_3)
  )
) +
  
  geom_col(
    width = 0.75
  ) +
  
  labs(
    title = "Symptom Scores — Time 1",
    x = "Sample",
    y = "Count of Symptom Scores",
    fill = "Symptom Score"
  ) +
  
  scale_fill_viridis_d(
    option = "D"
  ) +
  
  scale_x_discrete(
    labels = sample_labels
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(
      size = 10,
      vjust = 0.5,
      hjust = 1,
      angle = 90
    ),
    axis.title.y = element_text(
      size = 11
    ),
    axis.title.x = element_text(
      size = 11
    ),
    axis.text.y = element_text(
      size = 11
    ),
    plot.title = element_text(
      size = 11
    ),
    legend.text = element_text(
      size = 11
    ),
    legend.title = element_text(
      size = 11
    )
  )


# ============================================================
# Save graph
# ============================================================

ggsave(
  filename = file.path(
    output_dir,
    "Figure5-symptoms-t1.png"
  ),
  plot = p_symptoms,
  width = 8,
  height = 6,
  dpi = 300
)


# Display graph
print(p_symptoms)


# ============================================================
# Finished
# ============================================================

message(
  "Symptom T1 analysis complete."
)

message(
  "Dunn results saved to: ",
  file.path(
    output_dir,
    "Figure5-symptom-dunn-t1.csv"
  )
)

message(
  "Plot saved to: ",
  file.path(
    output_dir,
    "Figure5-symptoms-t1.png"
  )
)
