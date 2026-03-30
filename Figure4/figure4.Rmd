library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(tidyverse)
library(readxl)
library(pheatmap)
library(emmeans)
library(lme4)
library(ggpubr)


# Function for processing 96-well plates
# 1- Normalise each assay based on the mean blank well per plate + assay
# 2- Calculate mean per strain pseudorep (n = 4) per treatment per plate

process_assay <- function(data_dir,
                          plate_layout_file,
                          assay_meta_file,
                          experiment_name,
                          strain_order) {
  
  # ----------------------------
  # 1. Load metadata
  # ----------------------------
  
  plate_meta <- read_csv(plate_layout_file)
  
  assay_meta <- read_csv(assay_meta_file) %>%
    na.omit()
  
  # ----------------------------
  # 2. Read plate reader files
  # ----------------------------
  
  filenames <- list.files(
    path = data_dir,
    pattern = "\\.xlsx$",
    full.names = TRUE
  )
  
  # Extract and order plate numbers
  plate_numbers <- str_extract(basename(filenames), "(?<=plate)\\d+")
  order_idx <- order(as.integer(plate_numbers))
  
  filenames <- filenames[order_idx]
  plates    <- paste0("plate", plate_numbers[order_idx])
  
  # Read plates (8 × 12 from B42)
  plate_list <- setNames(
    lapply(
      filenames,
      read_excel,
      col_names = FALSE,
      range = anchored("B42", c(8, 12))
    ),
    plates
  )
  
  # ----------------------------
  # 3. Convert to long format
  # ----------------------------
  
  plate_df <- map2_dfr(
    plate_list,
    names(plate_list),
    ~ .x %>%
      as.data.frame() %>%
      setNames(1:12) %>%
      mutate(Row = LETTERS[1:8], .before = 1) %>%
      pivot_longer(
        cols = -Row,
        names_to = "Column",
        values_to = "OD600"
      ) %>%
      mutate(
        Column = as.integer(Column),
        Well   = paste0(Row, Column),
        plate  = .y
      )
  )
  
  # ----------------------------
  # 4. Merge metadata
  # ----------------------------
  
  merged_meta <- plate_df %>%
    left_join(plate_meta, by = "Well") %>%
    left_join(assay_meta, by = c("plate", "Assay_plate")) %>%
    mutate(
      Experiment = experiment_name,
      OD600 = as.numeric(OD600)
    )
  
  # ----------------------------
  # 5. Normalise to control (per plate)
  # ----------------------------
  
  controls <- merged_meta %>%
    filter(Strain == "Control") %>%
    group_by(plate, Assay_plate) %>%
    summarise(
      control_mean = mean(OD600, na.rm = TRUE),
      .groups = "drop"
    )
  
  normalised_meta <- merged_meta %>%
    left_join(controls, by = c("plate", "Assay_plate")) %>%
    mutate(OD600_norm = OD600 - control_mean)
  
  # ----------------------------
  # 6. Summarise per strain per plate
  # ----------------------------
  
  grouped_norm <- normalised_meta %>%
    filter(Strain != "Control") %>%   # remove control rows (optional but recommended)
    mutate(
      Strain = factor(Strain, levels = strain_order),
      plate = as.character(plate),
      Assay_plate = as.character(Assay_plate)
    ) %>%
    group_by(Strain, plate, Assay_plate, Mutation, Phenotype, Experiment) %>%
    summarise(
      mean_final = mean(OD600_norm, na.rm = TRUE),
      SD   = sd(OD600_norm, na.rm = TRUE),
      SE   = SD / sqrt(n()),
      n    = n(),
      .groups = "drop"
    )
  
  return(grouped_norm)
}

# ASSAY 1 ~~~~~~~~~~~~~~~~
grouped_norm_exp1 <- process_assay(
  data_dir = "data/Assay-1-151225",
  plate_layout_file = "data/PLATE_LAYOUT.csv",
  assay_meta_file = "data/assay_plate_metadata.csv",
  experiment_name = "Exp1-151225",
  strain_order = c(
    "W1413", "W1513", "W2013", "W2313",
    "W573", "W593",
    "W163_WT", "W163_GacA", "W163_GacS"
  )
)


# ASSAY 2 ~~~~~~~~~~~~~~~~
grouped_norm_exp2 <- process_assay(
  data_dir = "data/Assay-2-120126",
  plate_layout_file = "data/PLATE_LAYOUT_full-v2.csv",
  assay_meta_file = "data/assay_plate_metadata-exp2.csv",
  experiment_name = "Exp2-120126",
  strain_order = c(
    "W1413", "W1513", "W2013", "W2313",
    "W573", "W593",
    "W163_WT", "W163_GacA", "W163_GacS",
    "W163_GacS_GacS", "W163_GacS_EV"
  )
)

# ASSAY 3 ~~~~~~~~~~~~~~~~
grouped_norm_exp3 <- process_assay(
  data_dir = "data/Assay-3-190126",
  plate_layout_file = "data/PLATE_LAYOUT_full-v2.csv",
  assay_meta_file = "data/assay_plate_metadata-exp2.csv",
  experiment_name = "Exp3-190126",
  strain_order = c(
    "W1413", "W1513", "W2013", "W2313",
    "W573", "W593",
    "W163_WT",  "W163_GacA", "W163_GacS",
    "W163_GacS_GacS", "W163_GacS_EV"
  )
)

# ASSAY 4 ~~~~~~~~~~~~~~~~
grouped_norm_exp4 <- process_assay(
  data_dir = "data/Assay-4-260126",
  plate_layout_file = "data/PLATE_LAYOUT_full-v2.csv",
  assay_meta_file = "data/assay_plate_metadata-exp2.csv",
  experiment_name = "Exp4-260126",
  strain_order = c(
    "W1413", "W1513", "W2013", "W2313",
    "W573", "W593",
    "W163_WT",  "W163_GacA", "W163_GacS",
    "W163_GacA_GacA", "W163_GacA_EV"
  )
)

# Combine all experiments

Normalised_data_final <- bind_rows(
  grouped_norm_exp1,
  grouped_norm_exp2,
  grouped_norm_exp3,
  grouped_norm_exp4
)


# Generate boxplots for each phenotype, seperated by gacAS mutation status ~~~~~~~~~~~~~~~~~~~

filtered_df <- Normalised_data_final %>%
  filter(!is.na(Phenotype)) %>%
  filter(Strain != "W1413" | Phenotype != "M9_Glut_Gluc" | Experiment != "Exp3-190126") %>% # Removed as wells were not inoculated
  filter(Phenotype != "M9_Suc" | Experiment != "Exp1-151225") %>% # Removed because no growth in any well 
  filter(Strain != "W163_GacS_GacS",    # Removed as not enough replicates
         Strain != "W163_GacS_EV") %>%  # Removed as not enough replicates
  mutate(
    Genotype_type = case_when(
      Strain %in% c("W1513", "W2313", "W573", "W593", "W163_GacA", "W163_GacS") ~ "Mutant",
      Strain %in% c("W163_WT", "W1413", "W2013") ~ "Functional",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(log_mean = log(mean_final + 1))  # +1 if you have zeros

###############################################################
# Linear-mixed effects model  
###############################################################

# ----------------------------
# Prepare data
# ----------------------------
lmer_filtered_df <- filtered_df %>%
  filter(!is.na(Phenotype)) %>%
  mutate(
    Phenotype = factor(Phenotype),
    Genotype_type = factor(Genotype_type),
    Experiment = factor(Experiment),
    Strain = factor(Strain)
  )
write_csv(lmer_filtered_df, "figures/phenotype_data.csv")
# ----------------------------
# Mixed model
# ----------------------------
model_mixed <- lmer(
  mean_final ~ Phenotype * Genotype_type + (1 | Experiment),
  data = lmer_filtered_df
)

summary(model_mixed)

# ----------------------------
# Posthoc comparisons
# ----------------------------

# 1. Phenotype comparisons within genotype (not used for plot here, but retained)
emm <- emmeans(model_mixed, ~ Phenotype | Genotype_type)
functionalgroup_pairs <- as.data.frame(pairs(emm, adjust = "tukey"))
write_csv(functionalgroup_pairs, "figures/genotype_lmem.csv")
# 2. Functional vs Mutant within each phenotype (THIS is what we plot)
phenotype_emm <- emmeans(model_mixed, ~ Genotype_type | Phenotype)
write_csv(phenotype_pairs, "figures/phenotype_lmem.csv")
# ----------------------------------------------------------------
# Format original dataframe - subset to significant phenotypes only
# ----------------------------------------------------------------

sig_phenotypes <- plot_stats %>%
  filter(p.signif != "ns") %>%
  mutate(Phenotype = as.character(Phenotype)) %>%
  pull(Phenotype)

# Subset by significant results and enforce order
lmer_filtered_df <- lmer_filtered_df %>%
  filter(Phenotype %in% sig_phenotypes) %>%
  filter(!str_detect(Phenotype, "EToH")) %>% # Remove control media 
  mutate(Phenotype = factor(Phenotype, levels = c("M9_Glyc", "M9_Gluc", "M9_Ser_Gluc", "M9_Val_Gluc", "M9_Glut_Gluc", "M9_Pcoum_Gluc", "TSB_PEG"))
  )

# ----------------------------
# Format significance labels
# ----------------------------
plot_stats <- phenotype_pairs %>%
  filter(contrast == "Functional - Mutant") %>%
  separate(contrast, into = c("group1", "group2"), sep = " - ") %>%
  mutate(
    p.signif = case_when(
      p.value <= 0.0001 ~ "****",
      p.value <= 0.001  ~ "***",
      p.value <= 0.01   ~ "**",
      p.value <= 0.05   ~ "*",
      TRUE              ~ "ns"
    )
  )

# ----------------------------
# Add y positions per phenotype
# ----------------------------
plot_stats <- plot_stats %>%
  group_by(Phenotype) %>%
  mutate(
    y.position = max(
      lmer_filtered_df$mean_final[
        as.character(lmer_filtered_df$Phenotype) == as.character(Phenotype)
      ],
      na.rm = TRUE
    ) * 1.1
  ) %>%
  ungroup()

# ----------------------------
# Subset to signif stats only
# ----------------------------

plot_stats <- plot_stats %>%
  filter(Phenotype %in% sig_phenotypes) %>%
  filter(!str_detect(Phenotype, "EToH"))
  
# ----------------------------
# Facet labels and plot
# ----------------------------
facet_labels <- c(
  "M9_Glyc" = "M9 +\nGlycerol",
  "M9_Gluc" = "M9 +\nGlucose",
  "M9_Ser_Gluc" = "M9 +\nGlucose\n+ Serine",
  "M9_Glut_Gluc" = "M9 +\nGlucose\n+ Glutamate",
  "M9_Val_Gluc" = "M9 +\nGlucose\n+ Valine",
  "M9_Pcoum_Gluc" = "M9 +\nGlucose\n+ P-coumarin",
  "TSB_PEG" = "Water\npotential\nstress"
)


ggplot(lmer_filtered_df,
       aes(x = Genotype_type,
           y = mean_final)) +
  
  geom_boxplot(aes(fill = Genotype_type),
               outlier.shape = NA) +
  
  geom_jitter(
    aes(colour = Strain),
    width = 0.15,
    size = 1,
    alpha = 0.8
  ) +
  
  scale_fill_manual(values = c(
    "Functional" = "#FFD5D5",
    "Mutant"     = "#E3F4D7"
  )) +
  
  stat_pvalue_manual(
    plot_stats,
    label = "p.signif",
    tip.length = 0.01
  ) +
  
  facet_wrap(~ Phenotype,
             labeller = as_labeller(facet_labels),
             nrow = 1) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  
  theme_classic() +
  
  labs(
    x = NULL,
    y = "Normalised absorbance (OD600)",
    colour = "Strain"
  ) +
  
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    strip.text = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )


# ================================================================
# Generate heatmap for all exps: Plate-based growth assay analysis
# ================================================================

heatmap_df <- Normalised_data_final %>%
  filter(Strain != "W163_GacS_EV",
         Strain != "W163_GacS_GacS") %>%
  group_by(Strain, Phenotype) %>%
  summarise(
    mean_OD = mean(mean_final, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit() %>%
  filter(Phenotype != "M9_Frax",          # Removing exp with not enough replicates and controls
         Phenotype != "M9_Frax_Gluc",
         Phenotype != "M9_Pcoum",
         Phenotype != "M9_EToH100",
         Phenotype != "M9_EToH50",
         Phenotype != "M9" )

## Add swimming and biosurfactant
# Reload data
swim_df <- read_excel("data/phenotype_data_gacAS_240226.xlsx")
# Subset for swimming and exclude strains with Expression_vector
swim_df <- swim_df %>%
  filter(
    Assay == "Swimming",
    is.na(Expression_vector)
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Bio_rep = as.factor(Bio_rep)) %>%
  group_by(Strain_Genotype, Assay, experiment, Bio_rep) %>%
  summarise(
    mean_final = mean(Value),
    SD   = sd(Value),
    SE   = SD / sqrt(n()),
    .groups = "drop"
  ) %>%
  select(Strain_Genotype, Assay, mean_final) %>%
  magrittr::set_colnames(c("Strain", "Phenotype", "mean_OD"))

# Reload data
bs_df <- read_excel("data/phenotype_data_gacAS_240226.xlsx")
# Subset for swimming and exclude strains with Expression_vector
bs_df <- bs_df %>%
  filter(
    Assay == "Biosurfactant",
    is.na(Expression_vector)
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Bio_rep = as.factor(Bio_rep)) %>%
  group_by(Strain_Genotype, Assay, experiment, Bio_rep) %>%
  summarise(
    mean_final = mean(Value),
    SD   = sd(Value),
    SE   = SD / sqrt(n()),
    .groups = "drop"
  ) %>%
  select(Strain_Genotype, Assay, mean_final) %>%
  magrittr::set_colnames(c("Strain", "Phenotype", "mean_OD"))

all_df <- rbind(heatmap_df, swim_df, bs_df)

write_csv(all_df, "figures/heatmap_data.csv")

all_df_clean <- all_df %>%
  group_by(Strain, Phenotype) %>%
  summarise(
    mean_OD = mean(mean_OD, na.rm = TRUE),
    .groups = "drop"
  )

heatmap_mat <- all_df_clean %>%
  pivot_wider(
    names_from  = Phenotype,
    values_from = mean_OD
  ) %>%
  column_to_rownames("Strain") %>%
  as.matrix()

# Z-score by phenotype (column-wise scaling)
heatmap_mat_z <- scale(heatmap_mat)

pheatmap(
  heatmap_mat_z,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = NA,
)



