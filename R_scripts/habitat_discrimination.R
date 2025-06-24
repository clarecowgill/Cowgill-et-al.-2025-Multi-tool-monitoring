
# Loading necessary libraries
library(vegan)
library(tidyverse)
library(indicspecies)
library(pairwiseAdonis)

methods_week_pa <-read.csv('Data/all_week_pa.csv', row.names <-1)
week_ca_vegmet <-read.csv('Data/all_week_meta.csv', row.names <-2)

# Looks like this camera sample is a large outlier on NMDS though only includes Prunella modularis
# Scat samples only taken in confier and sample size too small anyway
week_met_no_outlier <- week_ca_vegmet[!grepl("^PO", rownames(week_ca_vegmet)), ]
week_met_no_outlier <-week_met_no_outlier[rownames(week_met_no_outlier) !="NCL08A_2023-06-03", ]
week_pa_no_outlier <-methods_week_pa[rownames(week_met_no_outlier), ]

# ______________________________________________________________________________

## Habitat distinction ----

pa_dat <-week_pa_no_outlier
meta_dat <-week_met_no_outlier

# Initialize empty lists for storing results
results_dispersion <-list()
results_permanova <-list()
results_nmds <-list()
results_indicators <-list()
stress_values <-data.frame(Method <-unique(meta_dat$type), Stress <-NA)

# Create an empty summary table
summary_table <-data.frame(
  Sample_Type <-character(0),
  PERMANOVA_R2 <-numeric(0),
  PERMANOVA_p_value <-numeric(0),
  Dispersion_Mean <-numeric(0),
  Dispersion_Var <-numeric(0),
  NMDS_Stress <-numeric(0),
  stringsAsFactors <-FALSE
)

# Dispersion analysis loop
for (method in unique(meta_dat$type)) {
  # Subset data for the current method
  method_samples <-meta_dat$type ==method
  method_data <-pa_dat[method_samples, ]
  habitat_groups <-meta_dat$habitat[method_samples]
  
  # Jaccard distance matrix
  method_dist <-vegdist(method_data, method <-"jaccard")
  
  # Homogeneity of dispersions
  method_var <-betadisper(method_dist, habitat_groups)
  
  # Store results
  results_dispersion[[method]] <-list(
    distances <-method_dist,
    dispersion <-method_var
  )
  
  # Print results
  print(paste("Dispersion analysis for sample type:", method))
  print(anova(method_var))
  print(permutest(method_var))
}

# PERMANOVA and pairwise PERMANOVA loop
for (method in unique(meta_dat$type)) {
  # Subset data for the current method
  method_samples <-meta_dat$type ==method
  method_data <-pa_dat[method_samples, ]
  habitat_groups <-meta_dat$habitat[method_samples]
  
  # Jaccard distance matrix
  method_dist <-vegdist(method_data, method <-"jaccard")
  
  # PERMANOVA
  method_adonis <-adonis2(method_dist ~ habitat_groups, 
                          data <-meta_dat[method_samples, ], 
                          weights <-as.vector(table(habitat_groups)),
                          permutations <-9999)
  
  # Store results
  results_permanova[[method]] <-method_adonis
  
  # Pairwise PERMANOVA
  pairwise_results <-pairwise.adonis(method_dist, habitat_groups)
  
  # Print results
  print(paste("PERMANOVA for sample type:", method))
  print(method_adonis)
  
  # Print pairwise PERMANOVA results
  print(paste("Pairwise PERMANOVA for sample type:", method))
  print(pairwise_results)
}


# NMDS analysis loop
for (method in unique(meta_dat$type)) {
  # Subset data for the current method
  method_samples <-meta_dat$type ==method
  method_data <-pa_dat[method_samples, ]
  habitat_groups <-meta_dat$habitat[method_samples]
  
  # Jaccard distance matrix
  method_dist <-vegdist(method_data, method <-"jaccard")
  
  # Perform NMDS
  method_nmds <-metaMDS(method_dist, dist <-"jaccard", k <-3, trymax <-1000)
  
  # Store the stress value in the stress_values data frame
  stress_values$Stress[stress_values$Method ==method] <-method_nmds$stress
  
  # Store NMDS results
  results_nmds[[method]] <-method_nmds
}

# Combine the NMDS results into one data frame for plotting
all_nmds_data <-do.call(rbind, lapply(names(results_nmds), function(method) {
  nmds <-results_nmds[[method]]
  method_samples <-meta_dat$type ==method
  data.frame(
    NMDS1 <-nmds$points[, 1],
    NMDS2 <-nmds$points[, 2],
    Method <-method,
    Habitat <-meta_dat$habitat[method_samples],
    location <-meta_dat$location[method_samples]
  )
}))

# Create the NMDS plot with stress annotation
ggplot(all_nmds_data, aes(x <-NMDS1, y <-NMDS2, color <-Habitat, shape <-location)) +
  geom_point(size <-3, alpha <-0.8) +
  stat_ellipse(aes(group <-Habitat), level <-0.95) +
  facet_wrap(~Method) +
  theme_minimal() +
  labs(title <-"NMDS Ordination by Sample Type",
       x <-"NMDS1", y <-"NMDS2") +
  geom_text(data <-stress_values, aes(x <--2, y <-2, label <-paste("Stress:", round(Stress, 3))),
            inherit.aes <-FALSE, size <-4, color <-"black", fontface <-"italic")

# Create NMDS plot with hulls
## Calculate convex hulls outside ggplot
hulls_data <-data.frame()

for (method in unique(all_nmds_data$Method)) {
  for (habitat in unique(all_nmds_data$Habitat)) {
    # Subset data for the current habitat and method
    subset_data <-all_nmds_data %>%
      filter(Method ==method, Habitat ==habitat)
    
    # Only calculate hull if there are more than 2 points
    if (nrow(subset_data) >=3) {
      # Calculate the convex hull for this habitat group
      chull_indices <-chull(subset_data$NMDS1, subset_data$NMDS2)
      chull_points <-subset_data[chull_indices, ]
      chull_points$group <-habitat  # Assign habitat as a group
      hulls_data <-bind_rows(hulls_data, chull_points)  # Append to hulls_data
    }
  }
}

## Create the NMDS plot with convex hulls
ggplot(all_nmds_data, aes(x <-NMDS1, y <-NMDS2, color <-Habitat)) +
  geom_point(size <-3, alpha <-0.8) +
  # Plot convex hulls using the pre-calculated hulls_data
  geom_polygon(data <-hulls_data, aes(x <-NMDS1, y <-NMDS2, fill <-Habitat, group <-Habitat), 
               alpha <-0.2, color <-"black") +
  facet_wrap(~Method, scales <-"fixed") +
  coord_fixed(ratio <-1) +
  theme(panel.background <-element_rect(fill <-'white'),
        axis.line.x <-element_line(colour <-'black', linewidth <-0.5, linetype <-'solid'),
        axis.line.y <-element_line(colour <-'black', linewidth <-0.5, linetype <-'solid'),
        axis.title.x <-element_text(margin <-unit(c(8, 0, 0, 0), "mm")),
        axis.title.y <-element_text(margin <-unit(c(0, 5, 0, 0), "mm")),
        axis.text.x <-element_text(colour <-"black"),
        axis.text.y <-element_text(colour <-"black"),
        plot.title <-element_text(face <-"bold", hjust <-0, colour <-"black"),
        plot.subtitle <-element_text(face <-"bold", hjust <-0, color <-"black", margin <-unit(c(2, 0, 0, 0), "mm")),
        text <-element_text(size <-14),
        legend.position <-"bottom",
        legend.direction <-"horizontal",
        legend.box <-"vertical",
        legend.key <-element_blank()) +
  labs(x <-"NMDS1", y <-"NMDS2") +
  geom_text(data <-stress_values, aes(x <--2, y <-2, label <-paste("Stress:", round(Stress, 3))),
            inherit.aes <-FALSE, size <-4, color <-"black", fontface <-"italic")

# Create a summary table with PERMANOVA, Dispersion, and NMDS Stress values
for (method in unique(meta_dat$type)) {
  permanova_result <-tryCatch(results_permanova[[method]], error <-function(e) NULL)
  dispersion_result <-tryCatch(results_dispersion[[method]], error <-function(e) NULL)
  nmds_result <-tryCatch(results_nmds[[method]], error <-function(e) NULL)
  
  # PERMANOVA
  permanova_R2 <-ifelse(!is.null(permanova_result) && !is.null(permanova_result$R2), permanova_result$R2[1], NA)
  permanova_p_value <-ifelse(!is.null(permanova_result) && !is.null(permanova_result$`Pr(>F)`), permanova_result$`Pr(>F)`[1], NA)
  
  # Dispersion
  dispersion_mean <-ifelse(!is.null(dispersion_result), mean(dispersion_result$dispersion$distances, na.rm <-TRUE), NA)
  dispersion_var <-ifelse(!is.null(dispersion_result), var(dispersion_result$dispersion$distances, na.rm <-TRUE), NA)
  
  # NMDS Stress
  nmds_stress <-ifelse(!is.null(nmds_result) && !is.null(nmds_result$stress), nmds_result$stress, NA)
  
  # Append to summary table
  summary_table <-rbind(summary_table, data.frame(
    Sample_Type <-method,
    PERMANOVA_R2 <-permanova_R2,
    PERMANOVA_p_value <-permanova_p_value,
    Dispersion_Mean <-dispersion_mean,
    Dispersion_Var <-dispersion_var,
    NMDS_Stress <-nmds_stress,
    stringsAsFactors <-FALSE
  ))
}

## View the summary table
print(summary_table)

# _____________________________________
### Additive PERMANOVA ----
# Understand if any method combos can detect better habitat distinction
# Load required libraries
library(vegan)
library(combinat) # For generating combinations

# Initialize results storage
all_combinations_results <- data.frame(
  Methods <-character(),
  R2 <-numeric(),
  p_value <-numeric(),
  stringsAsFactors <-FALSE
)

# Generate all possible combinations of methods
methods <- unique(meta_dat$type) # Extract method types
method_combinations <- lapply(seq_along(methods), function(i) {
  combn(methods, i, simplify <-FALSE)
}) %>% unlist(recursive <-FALSE)

# Iterate through each combination of methods
for (comb in method_combinations) {
  # Combine the names of methods into a single string for identification
  comb_name <- paste(comb, collapse <-" + ")
  
  # Select samples corresponding to the current combination of methods
  selected_samples <- meta_dat$type %in% comb
  combined_data <- pa_dat[selected_samples, ]
  habitat_groups <- meta_dat$habitat[selected_samples]
  
  # Ensure no missing values
  valid_rows <- complete.cases(habitat_groups, combined_data)
  combined_data <- combined_data[valid_rows, ]
  habitat_groups <- habitat_groups[valid_rows]
  
  # Skip if there are fewer than two habitats represented
  if (length(unique(habitat_groups)) < 2) {
    warning(paste("Combination", comb_name, "has fewer than 2 habitats. Skipping."))
    next
  }
  
  # Compute distance matrix
  dist_matrix <- vegdist(combined_data, method <-"jaccard")
  
  # Run PERMANOVA
  perm_result <- adonis2(dist_matrix ~ habitat_groups)
  
  # Store results
  all_combinations_results <- rbind(all_combinations_results, data.frame(
    Methods <-comb_name,
    R2 <-perm_result$R2[1],
    p_value <-perm_result$`Pr(>F)`[1],
    stringsAsFactors <-FALSE
  ))
}

sorted_results <- all_combinations_results[order(-all_combinations_results$R2), ]
print(head(sorted_results, 30)) # Top 15 combinations

# Create a bar plot for the top 15 combinations based on R2
ggplot(sorted_results[1:15, ], aes(x <-reorder(Methods, R2), y <-R2)) +
  geom_bar(stat <-"identity", fill <-"steelblue") +
  labs(x <-"Methods Combination", y <-"R²", title <-"Top 15 Methods Combinations by R²") +
  theme_minimal() +
  theme(axis.text.x <-element_text(angle <-45, hjust <-1)) +
  coord_flip() 


# Indicator species analysis loop
for (method in unique(meta_dat$type)) {
  # Subset data for the current method
  method_samples <-meta_dat$type ==method
  method_data <-pa_dat[method_samples, ]
  habitat_groups <-meta_dat$habitat[method_samples]
  
  # Indicator species analysis
  indicator_result <-multipatt(method_data, habitat_groups, func <-"r.g")
  
  # Store results
  results_indicators[[method]] <-indicator_result
  
  # Print summary
  print(paste("Indicator species analysis for method:", method))
  print(summary(indicator_result))
}


# _____________________________________
# Plot each NMDS seperately ----
methods <- unique(all_nmds_data$Method)

# Define a color palette for your habitats
habitat_colors <- c("B" <-"#31bf2c",
                    "C" <-"#439eba",
                    "F" <-"#e68d20",
                    "G" <-"#8332a8")

# Loop through methods and create separate plots
for (method in methods) {
  # Filter data for the current method
  method_data <- all_nmds_data %>% filter(Method ==method)
  
  # Separate G data from the other habitats
  non_G_data <- method_data %>% filter(Habitat !="G")
  G_data <- method_data %>% filter(Habitat =="G")
  
  # Get hull data for the method
  method_hulls <- hulls_data %>% filter(Method ==method)
  
  # Define stress for the current method
  stress_val <- stress_values %>% filter(Method ==method) %>% pull(Stress)
  
  # Create the plot
  p <- ggplot() +
    # Plot the non-G data first
    geom_point(data <-non_G_data, aes(x <-NMDS1, y <-NMDS2, color <-Habitat), size <-3, alpha <-0.8) +
    stat_ellipse(data <-non_G_data, aes(x <-NMDS1, y <-NMDS2, fill <-Habitat),
                 alpha <-0.3, geom <-"polygon", color <-NA, level <-0.95) +
  # Plot the G data last so it appears on top
    geom_point(data <-G_data, aes(x <-NMDS1, y <-NMDS2, color <-Habitat), size <-3, alpha <-0.7) +
    stat_ellipse(data <-G_data, aes(x <-NMDS1, y <-NMDS2, fill <-Habitat), alpha <-0.3, geom <-"polygon", color <-NA) +
    labs(x <-"NMDS1", y <-"NMDS2",
         subtitle <-paste("Stress:", round(stress_val, 3))) +  # Move stress to subtitle
    scale_color_manual(values <-habitat_colors) +
    scale_fill_manual(values <-habitat_colors) + 
    theme(panel.background <-element_rect(fill <-'white'),
          panel.border <-element_rect(colour <-"black", fill <-NA, size <-1),
          axis.line.x <-element_line(colour <-'black', linewidth <-0.5, linetype <-'solid'),
          axis.line.y <-element_line(colour <-'black', linewidth <-0.5, linetype <-'solid'),
          axis.title.x <-element_text(margin <-unit(c(8, 0, 0, 0), "mm")),
          axis.title.y <-element_text(margin <-unit(c(0, 5, 0, 0), "mm")),
          axis.text.x <-element_text(colour <-"black"),
          axis.text.y <-element_text(colour <-"black"),
          plot.title <-element_text(face <-"bold", hjust <-0.5, size <-14),
          plot.subtitle <-element_text(face <-"plain", hjust <-0, size <-10),  # Less bold, smaller, and left-aligned
          text <-element_text(size <-14),
          legend.position <-"bottom",
          legend.direction <-"horizontal",
          legend.box <-"vertical",
          legend.key <-element_blank(),
          aspect.ratio <-0.75) +
    coord_fixed(ratio <-1)
  
  # Print the plot
  print(p)
  
  # Save the plot to a file
  ggsave(filename <-paste0("NMDS_Method_jac_", method, ".png"), plot <-p, width <-5, height <-5)
}

