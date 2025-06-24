
library(vegan)
library(tidyverse)
library(fastDummies)
library(dplyr)
library(SYNCSA)
library(cluster)
library(dunn.test)
#_____________________________________________________________________________
# FD metrics ----
nodes <-read.csv('nodes.csv')
meta <-read.csv("all_week_meta.csv", row.names <-1)
pa <-read.csv('all_week_pa.csv', row.names <-1)

method_colors <- c(  "Aud" <-"#8de0da",
                     "Cam" <-"#ee881c",
                     "W"   <-"#435bf7",
                     "T"   <-"#4ab577",
                     "S"   <-"#f5e467"
)
## AvTD - average taxonomic distance ----

pa_filtered <- pa[rowSums(pa) > 2, ]

nodes_unique <- nodes %>%
  dplyr::select(species, genus, family, order, class, phylum) %>%
  distinct()

avtd_matrix <- nodes_unique %>%
  filter(species %in% colnames(pa)) %>%
  arrange(factor(species, levels <-colnames(pa)))

# Set species as row names
rownames(avtd_matrix) <- avtd_matrix$species
avtd_matrix <- avtd_matrix[, -1]

avtd_dist_matrix <- taxa2dist(as.matrix(avtd_matrix))
avtd_dist_matrix <- as.matrix(avtd_dist_matrix)

# Calculate Average Taxonomic Distinctiveness (ATD) for each sample
avtd_results <- apply(pa_filtered, 1, function(sample) {
  present_species <- which(sample == 1)
  
  if (length(present_species) > 1) {
    # Extract the distances for the present species
    distances <- avtd_dist_matrix[present_species, present_species]
    
    # Calculate the average taxonomic distance
    avtd <- mean(distances[upper.tri(distances)], na.rm <-TRUE)  # Use upper.tri to avoid redundancy
    return(avtd)
  } else {
    return(NA)  # Return NA if less than 2 species are present
  }
})

avtd_df <- data.frame(sample <-rownames(pa_filtered), AvTD <-avtd_results)
# Merge with metadata
avtd_df <- left_join(avtd_df, meta, by <-"sample")
# remove faecal samples as not enough for comparisons
avtd_df_filtered <- avtd_df %>%
  filter(type !="P")

avTD_boxplot <-ggplot(avtd_df_filtered, aes(x <-type, y <-AvTD, fill <-type)) +
  geom_boxplot(outlier.shape <-NA, color <-"black", alpha <-0.75) +  
  geom_jitter(width <-0.1, alpha <-0.5, size <-1) +
  scale_fill_manual(values <-method_colors) +
  labs(x <-NULL, y <-"AvTD") +
  theme_minimal() +
  theme(
    panel.border <-element_rect(color <-"black", fill <-NA, size <-1),
    panel.grid <-element_blank(),
    axis.title.y <-element_text(size <-10),
    legend.position <-"none"
  )

print(avTD_boxplot)

ggsave('avTD_boxplot.png', plot <-avTD_boxplot, height <-4, width <-4, dpi <-300, bg <-"white")

#_____________________________________________________________________________

## Trait dissimilarity matrix for FD and Rao's quadratic entropy ----

filtered_nodes <- nodes %>%
  filter(species %in% colnames(pa)) %>%
  arrange(factor(species, levels <-colnames(pa)))

# Step 2: Select relevant functional trait columns
trait_data <- filtered_nodes %>%
  dplyr::select(species, Mass, trophic_traits, Locomotion) %>% 
  mutate(Mass <-log(Mass + 1))
#mass_bin_levels <- c("XS", "S", "M", "L", "XL", "XXL")
#trait_matrix$Mass_bins <- factor(trait_matrix$Mass_bins, levels <-mass_bin_levels, ordered <-TRUE)

# Convert relevant columns to factors
# trait_data$Mass_bins <- factor(trait_matrix$Mass_bins, ordered <-TRUE)
trait_data$trophic_traits <- as.factor(trait_data$trophic_traits)
trait_data$Locomotion <- as.factor(trait_data$Locomotion)
# Check mass is numeric
trait_data$Mass <- as.numeric(trait_data$Mass)

## Rao's quadratic entropy ----
# Ensure species are row names in traits
traits <- trait_data %>%
  column_to_rownames(var <-"species")

# Step 2: Filter samples with 2 or fewer species
sample_species_counts <- rowSums(pa)  # Count species per sample
comm_filtered <- pa[sample_species_counts > 2, ]  # Keep samples with >2 species

# Step 3: Calculate RaoQ for filtered samples
rao_results <- rao.diversity(comm_filtered, traits, ord <-"metric", standardize <-TRUE)
rao_values <- rao_results$FunRao  # Extract RaoQ values

# Step 4: Convert rao_values to a data frame
rao_df <- data.frame(sample <-names(rao_values), RaoQ <-rao_values, row.names <-NULL)

# Step 6: Join with metadata
rao_df <- left_join(rao_df, meta, by <-"sample")

# Step 7: Remove faecal samples
rao_df <- rao_df %>%
  filter(type !="P")

# Step 8: Create the boxplot
rao_boxplot <- ggplot(rao_df, aes(x <-type, y <-RaoQ, fill <-type)) +
  geom_boxplot(outlier.shape <-NA, color <-"black", alpha <-0.75) +  
  geom_jitter(width <-0.1, alpha <-0.5, size <-1) +
  labs(x <-NULL, y <-"RaoQ") +
  scale_fill_manual(values <-method_colors) +
  theme_minimal() +
  theme(
    axis.title.y <-element_text(size <-10),
    panel.grid <-element_blank(),
    axis.line <-element_line(color <-"black"),
    panel.border <-element_rect(color <-"black", fill <-NA, size <-1),
    legend.position <-"none"
  )

# Step 9: Display the plot
print(rao_boxplot)

ggsave('rao_boxplot.png', plot <-rao_boxplot, height <-4, width <-4, dpi <-300, bg <-"white")

#_____________________________________________________________________________

fredun <- rao_results$FunRedundancy  # Extract Redundancy values

# Step 4: Convert rao_values to a data frame
fred_df <- data.frame(sample <-names(fredun), FRed <-fredun, row.names <-NULL)

# Step 5: Remove zero values
fred_df <- fred_df %>%
  filter(FRed > 0)

# Step 6: Join with metadata
fred_df <- left_join(fred_df, meta, by <-"sample")

# Step 7: Remove faecal samples
fred_df <- fred_df %>%
  filter(type !="P")

# Step 8: Create the boxplot
fred_boxplot <- ggplot(fred_df, aes(x <-type, y <-FRed, fill <-type)) +
  geom_boxplot(outlier.shape <-NA, color <-"black", alpha <-0.75) +  
  geom_jitter(width <-0.1, alpha <-0.5, size <-1) +
  labs(x <-NULL, y <-"FRed") +
  scale_fill_manual(values <-method_colors) +
  theme_minimal() +
  theme(
    axis.title.y <-element_text(size <-10),
    panel.grid <-element_blank(),
    axis.line <-element_line(color <-"black"),
    panel.border <-element_rect(color <-"black", fill <-NA, size <-1),
    legend.position <-"none"
  )

# Step 9: Display the plot
print(fred_boxplot)

ggsave('fred_boxplott.png', plot <-fred_boxplot, height <-4, width <-4, dpi <-300, bg <-"white")

#____________________________________________________________________________
# Normality checks

hist(avtd_df$AvTD) # right skewed
hist(rao_df$RaoQ) # still a little right-skewed after log transformation
hist(fred_df$FRed)

avtd_df <- avtd_df %>% filter(type !="P")

dunn.test(avtd_df$AvTD, avtd_df$type, method <-"bonferroni")
dunn.test(fred_df$FRed, fred_df$type, method <-"bonferroni")
dunn.test(rao_df$RaoQ, rao_df$type, method <-"bonferroni")

#______________________________________________________________________________
# PCoA ----

## with mass as a continuous vairbale ----
# Select the relevant variables from the nodes dataset
trait_data <- nodes %>% 
  dplyr::select(trophic_traits, Locomotion, Mass) %>% 
  mutate(Mass <-log(Mass + 1))

# Use dummy_cols to automatically create binary columns for all categorical variables
trait_data <- dummy_cols(trait_data)

# Drop the original categorical columns as they are now encoded as binary
trait_data <- trait_data %>%
  dplyr::select(-c(trophic_traits, Locomotion))  # Drop original categorical columns, we have the dummy variables now

# Assign species names as rownames
rownames(trait_data) <- nodes$species

# Subset and reorder matrices to ensure alignment of species in both matrices
species_in_both <- intersect(colnames(pa), rownames(trait_data))
pa_aligned <- pa[, species_in_both, drop <-FALSE]
trait_data_aligned <- trait_data[species_in_both, , drop <-FALSE]

# Scale the trait data (including continuous 'Mass')
trait_data_scaled <- scale(trait_data_aligned)

# --- Step 1: Data Preparation and Transformation ---

# Convert to a data frame to ensure Mass is treated correctly
trait_data_df <- as.data.frame(trait_data_scaled)

# --- Step 2: Calculate Gower Dissimilarity Matrix ---
# Calculate Gower dissimilarity matrix
gower_dist <- daisy(trait_data_scaled, metric <-"gower")

# --- Step 3: Perform PCoA ---
# PCoA on the Gower distance matrix
pcoa_result <- cmdscale(gower_dist, eig <-TRUE, k <-2)

# --- Step 4: Project Samples (Sites) ---
# Project samples using the presence/absence matrix
trait_by_sample <- as.matrix(pa_aligned) %*% as.matrix(pcoa_result$points)

# Convert to a data frame for plotting
pcoa_scores <- as.data.frame(trait_by_sample)
colnames(pcoa_scores) <- c("PC1", "PC2")
pcoa_scores$PC2 <- -pcoa_scores$PC2

# Merge with metadata for visualization
pcoa_scores <- cbind(pcoa_scores, meta[match(rownames(pcoa_scores), meta$sample), ])

# --- Step 5: Calculate Species Weights ---
# Calculate species positions as weighted averages of sample positions
species_weights <- t(as.matrix(pa_aligned)) %*% as.matrix(pcoa_scores[, c("PC1", "PC2")])
rownames(species_weights) <- colnames(pa_aligned)

# --- Step 6: Correlate Traits with PCoA Axes ---
# Correlate traits with PCoA axes
species_weights_scaled <- scale(species_weights)
trait_loadings <- cor(trait_data_scaled, species_weights_scaled)

# Prepare loadings for plotting
trait_loadings_df <- as.data.frame(trait_loadings)
trait_loadings_df$trait <- rownames(trait_loadings_df)

# --- Step 7: Scale Arrows for Visualization ---
# Scale arrows for visualization
arrow_scaling_factor <- max(abs(pcoa_scores$PC1), abs(pcoa_scores$PC2)) / 
  max(abs(trait_loadings_df$PC1), abs(trait_loadings_df$PC2))
trait_loadings_df$PC1_scaled <- trait_loadings_df$PC1 * arrow_scaling_factor
trait_loadings_df$PC2_scaled <- trait_loadings_df$PC2 * arrow_scaling_factor

# --- Step 8: Calculate Variance Explained by Each PCoA Axis ---
# Calculate the proportion of variance explained by each PCoA axis
variance_explained <- eigenvals(pcoa_result) / sum(eigenvals(pcoa_result))

# Create axis labels with % variance explained
pc_labels <- paste0("PC", 1:length(variance_explained), 
                    " (", round(variance_explained * 100, 2), "%)")

# --- Step 9: Add Custom Labels for Traits ---

# Custom labels for traits
custom_labels <- c(
  "trophic_traits_CV" <-"V",
  "trophic_traits_CI" <-"I",
  "trophic_traits_C" <-"C",
  "trophic_traits_O" <-"O",
  "trophic_traits_HF" <-"F",
  "trophic_traits_HG" <-"G",
  "trophic_traits_H" <-"H",
  "Locomotion_Fossorial" <-"Fos",
  "Locomotion_Aquatic" <-"Aqu",
  "Locomotion_Semi-aquatic" <-"SemAqu",
  "Locomotion_Semi-arboreal" <-"SemArb",
  "Locomotion_Arboreal" <-"Arb",
  "Locomotion_Terrestrial" <-"Ter",
  "Locomotion_Volant" <-"Aer",
  "Mass" <-"Mass"
)

# Add custom labels to the trait loadings dataframe
trait_loadings_df$custom_label <- custom_labels[match(rownames(trait_loadings_df), names(custom_labels))]

pcoa_scores <- pcoa_scores %>%
  filter(type !="P")

# Define custom color palette for types

# --- Step 10: Create PCoA Plot ---
pcoa_plot <- ggplot(pcoa_scores, aes(x <-PC1, y <-PC2, color <-type)) +
  geom_vline(xintercept <-0, color <-"gray", linetype <-"dashed", linewidth <-0.5) +
  geom_hline(yintercept <-0, color <-"gray", linetype <-"dashed", linewidth <-0.5) +
  geom_point(size <-3, alpha <-0.75) + 
  #stat_ellipse(aes(fill <-type), geom <-"polygon", alpha <-0.2, color <-NA) +
  geom_segment(data <-trait_loadings_df, aes(x <-0, y <-0, xend <-PC1_scaled, yend <-PC2_scaled),
               arrow <-arrow(type <-"closed", length <-unit(1, "mm")), linewidth <-0.2, color <-"black") +
  geom_text(data <-trait_loadings_df, aes(x <-PC1_scaled, y <-PC2_scaled, label <-custom_label),
            size <-5, color <-"black", vjust <--0.5) +
  labs(x <-pc_labels[1], y <-pc_labels[2], color <-"Sample Type") +
  theme_minimal() +
  scale_color_manual(name <-"Sample Type", values <-method_colors) +
  scale_fill_manual(values <-method_colors) +
  theme(
    legend.position <-"bottom",
    panel.grid.major <-element_blank(),
    panel.grid.minor <-element_blank(),
    #axis.line <-element_line(color <-"black"),
    panel.border <-element_rect(color <-"black", fill <-NA, size <-1)
  )

# Print the PCoA plot
print(pcoa_plot)

ggsave('PCoA_traits.png', plot <-pcoa_plot, height <-6.5, width <-8.5, dpi <-300, bg <-"white")

