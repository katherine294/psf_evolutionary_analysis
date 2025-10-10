# Load packages
library(ggplot2)
library(readr)
library(dplyr)
library(vcfR)
library(ape)

# Calculate pairwise distance for intra-site diversity
# Read in vcf file
vcf <- read.vcfR("~/W163a3b1_core_genome_masked.vcf")

# Convert vcf to binary data matrix 
geno <- extract.gt(vcf, element = "GT", as.numeric = TRUE) 
geno[is.na(geno)] <- 0  # Replace missing genotypes with 0

# Compute pairwsie SNP distance  - Hamming distance (counting SNP differences) to create a distance matrix.
snp_dist <- dist(t(geno), method = "manhattan")  # Manhattan distance counts differences
snp_dist
# Convert distance object into a matrix:
snp_dist_matrix <- as.matrix(snp_dist)
rownames(snp_dist_matrix) <- colnames(geno)
colnames(snp_dist_matrix) <- colnames(geno)


##### Intra-clade diversity #####
metadata <- read.csv("~/OneDrive - University of Birmingham/Genotype-phenotype paper/Data for figures/metadata_clade.csv")
metadata <- metadata[c(1,2)]
metadata <- na.omit(metadata)

# Step 2: Define the function to compute mean intra-group SNP distance
calculate_intra_clade_diversity <- function(group_id, snp_dist_matrix, metadata) {
  isolates <- metadata$Strain[metadata$Clade == group_id]
  
  # Keep only strains that exist in the SNP distance matrix
  isolates <- isolates[isolates %in% rownames(snp_dist_matrix)]
  
  if (length(isolates) < 2) {
    return(NA)  # Return 0 instead of NA for singletons
  }
  
  sub_matrix <- snp_dist_matrix[isolates, isolates, drop = FALSE]
  mean_distance <- mean(sub_matrix[upper.tri(sub_matrix)])
  
  return(mean_distance)
}

# Step 3: Apply to each site_tree_tissue group
intra_clade_diversity <- metadata %>%
  distinct(Clade) %>%
  rowwise() %>%
  mutate(Mean_SNP_Divergence = calculate_intra_clade_diversity(Clade, snp_dist_matrix, metadata)) %>%
  ungroup()

# Calculate summary statistics
within_clade_summary <- data.frame(
  Mean = mean(intra_clade_diversity$Mean_SNP_Divergence),
  SD = sd(intra_clade_diversity$Mean_SNP_Divergence),
  N_Pairs = length(intra_clade_diversity$Mean_SNP_Divergence)
)

#### Divergence BETWEEN clades #####

# Get list of unique clades
clades <- unique(metadata$Clade)

# Initialize empty dataframe
clade_divergence <- data.frame()

# Loop through all pairwise clade combinations
for (i in 1:(length(clades)-1)) {
  for (j in (i+1):length(clades)) {
    clade1 <- clades[i]
    clade2 <- clades[j]
    
    isolates1 <- metadata$Strain[metadata$Clade == clade1]
    isolates2 <- metadata$Strain[metadata$Clade == clade2]
    
    # Keep only strains that exist in the distance matrix
    isolates1 <- isolates1[isolates1 %in% rownames(snp_dist_matrix)]
    isolates2 <- isolates2[isolates2 %in% rownames(snp_dist_matrix)]
    
    # Subset the SNP matrix
    sub_matrix <- snp_dist_matrix[isolates1, isolates2, drop=FALSE]
    
    # Compute mean divergence
    mean_div <- mean(sub_matrix, na.rm=TRUE)
    
    # Store result
    clade_divergence <- rbind(clade_divergence, data.frame(
      Clade1 = clade1,
      Clade2 = clade2,
      Mean_SNP_Divergence = mean_div
    ))
  }
}

# View results
print(clade_divergence)

# Rename and reshape intra-clade data
intra_formatted <- intra_clade_diversity %>%
  mutate(Clade1 = Clade, Clade2 = Clade) %>%
  select(Clade1, Clade2, Mean_SNP_Divergence)


# Calculate summary statistics
between_clade_summary <- data.frame(
  Mean = mean(clade_divergence$Mean_SNP_Divergence),
  SD = sd(clade_divergence$Mean_SNP_Divergence),
  N_Pairs = length(clade_divergence$Mean_SNP_Divergence)
)

# Combine with clade divergence
full_clade_divergence <- bind_rows(clade_divergence, intra_formatted)


subset(full_clade_divergence, (Clade1 == "6" & Clade2 == "7") | (Clade1 == "7" & Clade2 == "6"))

#### Adding 6 and 7 reverse
# Extract the row with (7,6)
row_7_6 <- subset(full_clade_divergence, Clade1 == "7" & Clade2 == "6")
# Create the reverse pair (6,7) with same value
row_6_7 <- row_7_6
row_6_7$Clade1 <- "6"
row_6_7$Clade2 <- "7"
row_6_7$Clade1 <- as.numeric(row_6_7$Clade1)
row_6_7$Clade2 <- as.numeric(row_6_7$Clade2)
# Add to the dataframe
full_clade_divergence <- bind_rows(full_clade_divergence, row_6_7)

# Order clades
levels_order <- as.character(1:10)

full_clade_divergence$Clade1 <- factor(full_clade_divergence$Clade1, levels = levels_order)
full_clade_divergence$Clade2 <- factor(full_clade_divergence$Clade2, levels = levels_order)

# Add numeric versions of factors for filtering
full_clade_divergence <- full_clade_divergence %>%
  mutate(
    Clade1_num = as.numeric(as.character(Clade1)),
    Clade2_num = as.numeric(as.character(Clade2))
  )

# Keep only diagonal and upper triangle (for example)
plot_data <- full_clade_divergence %>%
  filter(Clade2_num >= Clade1_num)  # Keep diagonal and upper triangle only

heatmap <- ggplot(plot_data, aes(x = Clade1, y = Clade2, fill = Mean_SNP_Divergence)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "plasma", na.value = "grey90") +
  #scale_fill_gradient(low = "white", high = "steelblue", na.value = "grey90") +
  theme_minimal() +
  coord_fixed() + 
  geom_text(aes(label = round(Mean_SNP_Divergence, 1)), size = 3) +
  labs(title = "Within- and Between-Clade SNP Diversity",
       x = "Clade 1", y = "Clade 2", fill = "Mean SNP Divergence") +
  scale_y_discrete(limits = rev(levels(full_clade_divergence$Clade2))) 

heatmap
