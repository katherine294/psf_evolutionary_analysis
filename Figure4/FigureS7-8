############################################################
# Supplementary 7 and 8 
# Swarming Motility Analysis – GacA/S Mutants and WT Strains
# Author: Katherine
############################################################

## ──────────────────────────────────────────────────────────
## Libraries
## ──────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)
library(multcompView)
library(stringr)


############################################################
## MUTANT STRAINS
############################################################

# Load data
df <- read_excel("data/phenotype_data_gacAS_240226.xlsx")

# Subset for swarming assay (W163a3b1 only)
df <- df %>%
  filter(
    Assay == "Swarming",
    Strain == "W163a3b1"
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Bio_rep = as.factor(Bio_rep),
    Strain_mutant = factor(
      Strain_mutant,
      levels = c(
        "W163a3b1_WT",
        "W163a3b1_WT_EV",
        "W163a3b1_WT_pBBR1_GacA",
        "W163a3b1_WT_pBBR1_GacS",
        "W163a3b1_GacA",
        "W163a3b1_GacA_EV",
        "W163a3b1_GacA_pBBR1_GacA",
        "W163a3b1_GacS",
        "W163a3b1_GacS_EV",
        "W163a3b1_GacS_pBBR1_GacS"
      )
    ),
    Genotype = factor(Genotype, levels = c("WT", "GacA", "GacS"))
  )

# ANOVA + Tukey
anova_model <- aov(Value ~ Strain_mutant, data = df)
summary(anova_model)

tukey <- TukeyHSD(anova_model)
print(tukey)

# Compact Letter Display
cld <- multcompLetters4(anova_model, tukey)

# Summary statistics for CLD placement
Tk <- df %>%
  group_by(Strain_mutant) %>%
  summarise(
    mean  = mean(Value, na.rm = TRUE),
    quant = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

Tk$cld <- as.factor(cld$Strain_mutant$Letters)

# Merge genotype for faceting
Tk <- Tk %>%
  left_join(
    df %>% select(Strain_mutant, Genotype),
    by = "Strain_mutant"
  )

# Facet labels
facet_labels <- c(
  "WT"   = "WT",
  "GacA" = "gacAΔ",
  "GacS" = "gacSΔ"
)

# Plot
ggplot(df, aes(x = Strain_mutant, y = Value)) +
  geom_boxplot() +
  geom_jitter(
    aes(colour = Bio_rep),
    width = 0.1,
    size = 1.5,
    alpha = 0.9
  ) +
  geom_text(
    data = Tk,
    aes(y = quant + 13, label = cld),
    size = 4
  ) +
  facet_wrap(
    ~ Genotype,
    scales = "free_x",
    labeller = as_labeller(facet_labels)
  ) +
  scale_x_discrete(
    labels = c(
      "",
      "pBBR1 EV",
      "pBBR1 + GacA",
      "pBBR1 + GacS",
      "gacAΔ",
      "gacAΔ + EV",
      "gacAΔ + GacA",
      "gacSΔ",
      "gacSΔ + EV",
      "gacSΔ + GacS"
    )
  ) +
  labs(
    x = "Expression vector",
    y = "Diameter (mm)",
    colour = str_wrap("Biological Replicate", width = 10),
    title = "Swarming motility of W163a3b1 Psf strains"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title   = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(size = 10),
    strip.text   = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  ylim(0, max(df$Value, na.rm = TRUE) + 13)



############################################################
## WILDTYPE STRAINS
############################################################

# Reload data
df <- read_excel("data/phenotype_data_gacAS_240226.xlsx")

# Subset for swarming (exclude mutants)
df <- df %>%
  filter(
    Assay == "Swarming",
    !Strain_Genotype %in% c(
      "W163a3b1_GacA",
      "W163a3b1_GacS"
    )
  ) %>%
  mutate(
    Value = as.numeric(Value),
    Bio_rep = as.factor(Bio_rep),
    Strain_mutant = factor(
      Strain_mutant,
      levels = c(
        "W163a3b1_WT",
        "W163a3b1_WT_EV",
        "W163a3b1_WT_pBBR1_GacA",
        "W163a3b1_WT_pBBR1_GacS",
        "W2313a1b3_WT",
        "W2313a1b3_WT_EV",
        "W2313a1b3_WT_pBBR1_GacA",
        "W1513a1b2_WT",
        "W1513a1b2_WT_EV",
        "W1513a1b2_WT_pBBR1_GacS",
        "W573a5b1_WT",
        "W593a5b1_WT"
      )
    ),
    Strain = factor(
      Strain,
      levels = c(
        "W163a3b1",
        "W2313a1b3",
        "W1513a1b2",
        "W573a5b1",
        "W593a5b1"
      )
    ),
    Strain_mutant = fct_drop(Strain_mutant)
  )

# ANOVA + Tukey
anova_model <- aov(Value ~ Strain_mutant, data = df)
summary(anova_model)

tukey <- TukeyHSD(anova_model)
print(tukey)

# Compact Letter Display
cld <- multcompLetters4(anova_model, tukey)

# Summary statistics
Tk <- df %>%
  group_by(Strain_mutant) %>%
  summarise(
    mean  = mean(Value, na.rm = TRUE),
    quant = quantile(Value, 0.75, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean))

Tk$cld <- as.factor(cld$Strain_mutant$Letters)

Tk <- Tk %>%
  left_join(
    df %>% select(Strain_mutant, Strain),
    by = "Strain_mutant"
  ) %>%
  distinct()

# Facet labels
facet_labels <- c(
  "W163a3b1"  = "W163 WT",
  "W2313a1b3" = "W2313 WT gacA",
  "W1513a1b2" = "W1513 WT gacS",
  "W573a5b1"  = "W573 WT gacS",
  "W593a5b1"  = "W593 WT gacA"
)

# Plot
ggplot(df, aes(x = Strain_mutant, y = Value)) +
  geom_boxplot() +
  geom_jitter(
    aes(colour = Bio_rep),
    width = 0.1,
    size = 1.5,
    alpha = 0.9
  ) +
  geom_text(
    data = Tk,
    aes(y = quant + 4, label = cld),
    size = 4
  ) +
  facet_wrap(
    ~ Strain,
    nrow = 1,
    scales = "free_x",
    labeller = as_labeller(facet_labels),
    drop = TRUE
  ) +
  scale_x_discrete(
    labels = c(
      "None",
      "pBBR1 EV",
      "pBBR1 + GacA",
      "pBBR1 + GacS",
      "gacAΔ",
      "gacAΔ + EV",
      "gacAΔ + GacA",
      "gacSΔ",
      "gacSΔ + EV",
      "gacSΔ + GacS",
      "",
      ""
    )
  ) +
  labs(
    x = "Expression vector",
    y = "Halo (mm)",
    colour = str_wrap("Biological Replicate", width = 10),
    title = "Swarming motility of wild-type Psf strains"
  ) +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title   = element_text(size = 10),
    axis.text.y  = element_text(size = 10),
    plot.title   = element_text(size = 10),
    strip.text   = element_text(size = 10),
    legend.text  = element_text(size = 10),
    legend.title = element_text(size = 10)
  ) +
  ylim(0, max(df$Value, na.rm = TRUE) + 3)
