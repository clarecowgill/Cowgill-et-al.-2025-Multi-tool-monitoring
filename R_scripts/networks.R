# Network analyses ----

library(vegan)
library(igraph)
library(dplyr)
library(cheddar)
library(NetIndices)
library(ggplot2)
library(cowplot)
library(scales)

#______________________________________________________________________________
# Get overall trophic levels ----

## Data sorting ----
nodes <- read.csv('nodes.csv')

nodes$node <- gsub(" ", "_", nodes$node)

node_names <- as.character(nodes$node)
trophic_links <- read.csv('trophic_links.csv')
# Remove duplicate links
trophic_links <- unique(trophic_links)
# remove dodgey looking links (e.g. chaffinch eating sparrowhawk)
trophic_links <- trophic_links %>%
  filter(consumer != "Fringilla coelebs") %>%
  unique()
# Define properties
properties <- list(title <- "NCL database")

trophic_links$consumer <- gsub(" ", "_", trophic_links$consumer)
trophic_links$resource <- gsub(" ", "_", trophic_links$resource)

trophic_links <- trophic_links %>%
  filter(consumer %in% nodes$node & resource %in% nodes$node)%>%
  distinct(consumer, resource, .keep_all <- TRUE)

#_______________
### Check for whether there are any missing nodes from trophic links ----

# Get all species that appear as consumers or resources in the trophic links
trophic_species <- unique(c(trophic_links$consumer, trophic_links$resource))

# Find species from the node list that are not in the trophic links
missing_species <- setdiff(nodes$node, trophic_species)

# Print missing species
if (length(missing_species) > 0) {
  cat("Species not present in trophic links:\n")
  print(missing_species)
} else {
  cat("All species in the node list have trophic links.")
}
#_______________

## Create the community ----
NCL.community <- Community(nodes, properties, trophic_links)

# Examining the different aspects of the community
NumberOfNodes(NCL.community)
NumberOfTrophicLinks(NCL.community)
# Node properties
NPS(NCL.community)
# Trophic links
TLPS(NCL.community)

# Find PreyAveragedTL
link_pos <- as.data.frame(TrophicLevels(NCL.community))
link_pos$rownames <- rownames(link_pos)
link_pos <- link_pos['PreyAveragedTL']
# Add 1 to all PreyAveragedTL values to shift the levels as the food web does not contain producers
link_pos$PreyAveragedTL <- link_pos$PreyAveragedTL + 1
link_pos$species <- rownames(link_pos)
mean(link_pos$PreyAveragedTL)
#______________________________________________________________________________
# Comparing network metrics by sample ----

## Data sorting ----
ncl_pa <- read.csv('cooccurrence_pa.csv', row.names <- 1)
colnames(ncl_pa)[colnames(ncl_pa) == "Apodemus"] <- "Apodemus_sylvaticus"

ncl_meta <- read.csv('cooccurrence_meta.csv')
nodes$node <- gsub(" ", "_", nodes$node)

trophic_links <- read.csv('trophic_links.csv')

# Remove duplicate links
trophic_links <- unique(trophic_links)
# remove dodgey looking links (e.g. chaffinch eating sparrowhawk)
trophic_links <- trophic_links %>%
  filter(consumer != "Fringilla coelebs") %>%
  unique()
# for consistent naming
trophic_links$consumer <- gsub(" ", "_", trophic_links$consumer)
trophic_links$resource <- gsub(" ", "_", trophic_links$resource)
trophic_links <- trophic_links %>%
  filter(consumer %in% nodes$node & resource %in% nodes$node)%>%
  distinct(consumer, resource, .keep_all <- TRUE)
# Remove self-links
trophic_links <- trophic_links %>%
  filter(consumer != resource)

write.csv(trophic_links, 'trophic_links_filtered.csv')
#_______________

## Calculate metrics ----

# Function to calculate PreyAveragedLength for detected species
calculate_chain_averaged_length <- function(graph) {
  if (gsize(graph) == 0) return(NA)  # If no links, return NA

  # Get the distances matrix
  distance_matrix <- distances(graph)

  # Set the distance for basal species to 1
  # Assuming basal species are those with degree 0 (no incoming edges)
  basal_species <- V(graph)[degree(graph, mode <- "in") == 0]
  distance_matrix[basal_species, basal_species] <- 1

  # Calculate average chain length (excluding infinite distances)
  avg_chain_length <- mean(distance_matrix[distance_matrix != Inf], na.rm <- TRUE)

  return(avg_chain_length)
}

# Function to calculate robustness with random removals
calculate_robustness <- function(graph, num_iterations <- 100) {
  if (gorder(graph) == 0) return(NA)  # If no nodes, return NA

  species_remaining <- gorder(graph)
  robustness_values <- numeric(num_iterations)

  for (iteration in 1:num_iterations) {
    g_copy <- graph
    removal_order <- sample(V(g_copy)$name)  # Random

    for (i in seq_along(removal_order)) {
      node_to_remove <- which(V(g_copy)$name == removal_order[i])
      g_copy <- delete_vertices(g_copy, node_to_remove)

      # Identify isolated nodes
      isolated_nodes <- V(g_copy)[degree(g_copy) == 0]
      if (length(isolated_nodes) > 0) {
        g_copy <- delete_vertices(g_copy, isolated_nodes)
      }

      if (gorder(g_copy) / species_remaining <<- 0.5) {
        robustness_values[iteration] <- i / species_remaining
        break
      }
    }
  }

  return(mean(robustness_values, na.rm <- TRUE))
}

# Function to calculate PreyAveragedLength (average chain length)
calculate_chain_averaged_length <- function(graph) {
  if (gsize(graph) == 0) return(NA)

  distance_matrix <- distances(graph)

  # Set distance for basal species to 1 (optional heuristic)
  basal_species <- V(graph)[degree(graph, mode <- "in") == 0]
  distance_matrix[basal_species, basal_species] <- 1

  avg_chain_length <- mean(distance_matrix[distance_matrix != Inf], na.rm <- TRUE)

  return(avg_chain_length)
}

network_metrics <- data.frame(sample <- character(),
                              type <- character(),
                              detected_species <- character(),  # Change to character
                              num_links <- numeric(),
                              linkage_density <- numeric(),
                              connectance <- numeric(),
                              nestedness <- numeric(),
                              modularity <- numeric(),
                              coherence <- numeric(),
                              robustness <- numeric(),
                              average_chain_length <- numeric(),
                              mean_trophic_level <- numeric(),
                              stringsAsFactors <- FALSE)

# Loop through each sample
for (s in rownames(ncl_pa)) {
  detected_species <- colnames(ncl_pa)[ncl_pa[s, ] == 1]
  
  # Skip samples with no detected species
  if (length(detected_species) == 0) {
    warning(paste("No detected species in sample", s, "- skipping."))
    next
  }
  
  # Debug: Print detected species
  print(paste("Detected species in sample", s, ":", paste(detected_species, collapse <- ", ")))
  
  # Get trophic levels
  detected_trophic_levels <- link_pos$PreyAveragedTL[link_pos$species %in% detected_species]
  
  # Debug: Print trophic levels
  print(paste("Trophic levels for detected species:", paste(detected_trophic_levels, collapse <- ", ")))
  
  # Remove NAs and compute mean
  trophic_levels <- na.omit(detected_trophic_levels)
  mean_trophic_level <- if (length(trophic_levels) == 0) NA else mean(trophic_levels)
  
  # Filter links
  sample_links <- trophic_links %>%
    filter(resource %in% detected_species & consumer %in% detected_species)
  
  # Skip if there are no links
  if (nrow(sample_links) == 0) {
    warning(paste("No valid links in sample", s, "- skipping."))
    next
  }
  
  # Create vertices data frame
  vertices_df <- data.frame(name <- gsub(" ", "_", detected_species), stringsAsFactors <- FALSE)
  
  # Build graph
  g <- graph_from_data_frame(sample_links, directed <- TRUE, vertices <- vertices_df)
  
  # Community object
  community <- Community(
    nodes <- data.frame(node <- V(g)$name),
    trophic.links <- sample_links,
    properties <- list(title <- s)
  )
  
  # Metrics
  num_links <- gsize(g)
  linkage_density <- ifelse(gorder(g) > 0, mean(degree(g)), NA)
  
  possible_links <- gorder(g) * (gorder(g) - 1)
  connectance <- ifelse(possible_links > 0, gsize(g) / possible_links, NA)
  
  adj_matrix <- as.matrix(as_adjacency_matrix(g))
  nestedness <- if (nrow(adj_matrix) > 1 && ncol(adj_matrix) > 1) nestednodf(adj_matrix)$statistic["NODF"] else NA
  
  modularity_value <- if (gsize(g) > 0) modularity(cluster_walktrap(g)) else NA
  
  sp_matrix <- distances(g)
  coherence_value <- if (num_links > 0) var(sp_matrix[sp_matrix != Inf], na.rm <- TRUE) else NA
  
  robustness_value <- calculate_robustness(g)
  avg_chain_length <- calculate_chain_averaged_length(g)
  
  # Look up sample type safely
  sample_type <- ncl_meta$type[ncl_meta$sample == s]
  if (length(sample_type) == 0) sample_type <- NA
  
  # Mismatch warning
  detected_count <- length(detected_species)
  expected_count <- sum(ncl_pa[s, ])
  if (detected_count != expected_count) {
    warning(paste("Mismatch in sample:", s, "- Detected species:", detected_count,
                  "Expected species:", expected_count))
  }
  
  # Assemble result row
  new_row <- data.frame(
    sample <- s,
    type <- sample_type,
    detected_species <- paste(detected_species, collapse <- ", "),
    num_links <- num_links,
    linkage_density <- linkage_density,
    connectance <- connectance,
    nestedness <- nestedness,
    modularity <- modularity_value,
    coherence <- coherence_value,
    robustness <- robustness_value,
    average_chain_length <- avg_chain_length,
    mean_trophic_level <- mean_trophic_level,
    stringsAsFactors <- FALSE
  )
  
  # Append
  network_metrics <- rbind(network_metrics, new_row)
}


#_______________
### Proportion of samples with no trophic links ----
# Count samples with no links for each type
no_links_summary <- data.frame(type <- character(),
                               total_samples <- integer(),
                               no_links_count <- integer(),
                               proportion_no_links <- numeric(),
                               stringsAsFactors <- FALSE)

for (sample_type in unique(ncl_meta$type)) {
  total_samples <- sum(ncl_meta$type == sample_type)

  # Count samples with no links
  no_links_count <- sum(network_metrics$num_links[network_metrics$type == sample_type] == 0)

  proportion_no_links <- no_links_count / total_samples

  # Append to summary dataframe
  no_links_summary <- rbind(no_links_summary, data.frame(type <- sample_type,
                                                         total_samples <- total_samples,
                                                         no_links_count <- no_links_count,
                                                         proportion_no_links <- proportion_no_links))
}

# Print summary
print(no_links_summary)

#_______________
## Plot boxplots of metrics ----
filtered_network_metrics <- network_metrics
filtered_network_metrics <- network_metrics[network_metrics$num_links > 0, ]
filtered_network_metrics <- filtered_network_metrics %>%
  filter(type != "P")

library(dunn.test)

dunn.test(filtered_network_metrics$linkage_density, filtered_network_metrics$type, method <- "bonferroni")
dunn.test(filtered_network_metrics$connectance, filtered_network_metrics$type, method <- "bonferroni")
dunn.test(filtered_network_metrics$robustness, filtered_network_metrics$type, method <- "bonferroni")
dunn.test(filtered_network_metrics$mean_trophic_level, filtered_network_metrics$type, method <- "bonferroni")


metrics <- c("linkage_density", "connectance","robustness",
             "modularity", "average_chain_length", "mean_trophic_level")

method_colors <- c(  "Aud" <- "#8de0da",
                     "Cam" <- "#ee881c",
                     "W"   <- "#435bf7",
                     "T"   <- "#4ab577",
                     "S"   <- "#f5e467"
                  )

plot_titles <- c("D)", "E)", "F)", "G)", "H", "I")

axes_titles <- c("Linkage Density", "Connectance", "Robustness (R50)", "Modularity", "Av Chain Length", "Mean Trophic Level")

plot_list <- list()

for (i in seq_along(metrics)) {
  p <- ggplot(filtered_network_metrics, aes_string(x <- "type", y <- metrics[i], fill <- "type")) +
    geom_boxplot(width <- 0.8, outlier.shape <- NA, alpha <- 0.8) +
    geom_jitter(width <- 0.1, alpha <- 0.5, size <- 1) +
    scale_fill_manual(values <- method_colors) +
    labs(title <- plot_titles[i], x <- NULL, y <- axes_titles[i]) +
    theme_minimal() +
    theme(
      panel.border <- element_rect(color <- "black", fill <- NA, size <- 1),
      panel.grid.major <- element_blank(),
      panel.grid.minor <- element_blank(),
      axis.title.y <- element_text(size <- 16),
      legend.position <- "none"
      )

  plot_list[[i]] <- p
}

combined_plot <- plot_grid(plotlist <- plot_list, ncol <- 2)
print(combined_plot)
ggsave('metric_plots.png', combined_plot, dpi<- 300, bg <- 'white')


#____________________________________________________________________________
# Overlaying all networks ----
## Define method-specific colors
method_colors_aud <- c(
  "Aud" <- "#8de0da",
  "Cam" <- "#ededed",
  "P"   <- "#ededed",
  "S"   <- "#ededed",
  "T"   <- "#ededed",
  "W"   <- "#ededed"
)

method_colors_cam <- c(
  "Aud" <- "#ededed",
  "Cam" <- "#ee881c",
  "P"   <- "#ededed",
  "S"   <- "#ededed",
  "T"   <- "#ededed",
  "W"   <- "#ededed"
)

method_colors_dna <- c(
  "Aud" <- "#ededed",
  "Cam" <- "#ededed",
  "P"   <- "#ededed",
  "S"   <- "#f5c905",
  "T"   <- "#2cd486",
  "W"   <- "#435bf7"
)

method_colors_all <- c(
  "Aud" <- "#8de0da",
  "Cam" <- "#ee881c",
  "P"   <- "#ededed",
  "S"   <- "#f5c905",
  "T"   <- "#2cd486",
  "W"   <- "#435bf7"
)

method_colors_all_grey <- c(
  "Aud" <- "grey",
  "Cam" <- "grey",
  "P"   <- "grey",
  "S"   <- "grey",
  "T"   <- "grey",
  "W"   <- "grey"
)

method_colors <- method_colors_aud

## Set positions -----------------------------------------------------------
detected_species <- colnames(ncl_pa)
set.seed(42)

# Get trophic levels and initialize positions
vertices_df <- data.frame(name <- detected_species) %>%
  left_join(link_pos, by <- c("name" <- "species"))

fixed_x_positions <- numeric(length(detected_species))
names(fixed_x_positions) <- detected_species

# Automated x-position assignment
for (tl in unique(vertices_df$PreyAveragedTL)) {
  species_at_tl <- vertices_df$name[vertices_df$PreyAveragedTL == tl]
  available_positions <- seq(-4.5, 4.5, by <- 0.1)
  sampled_positions <- sample(available_positions, length(species_at_tl), replace <- FALSE)
  fixed_x_positions[species_at_tl] <- sampled_positions
}

# Manual position overrides
manual_positions <- c(
  "Corvus_corax" <- -4, 
  "Corvus_corone" <- -3, 
  "Strix_aluco" <- -2.6, 
  "Numenius_arquata" <- -3.4, 
  "Garrulus_glandarius" <- -1.8, 
  "Buteo_buteo" <- -0.6, 
  "Accipiter_gentilis" <- 0.2, 
  "Vulpes_vulpes" <- 3.8, 
  "Dendrocopos_major" <- 0, 
  "Martes_martes" <- 1.2, 
  "Meles_meles" <- 2,
  "Sciurus_vulgaris" <- 3.2,
  "Phasianidae" <- 0.6,
  
  "Parus_major" <- 1,
  "Rattus_norvegicus" <- -0.8,
  "Turdus_philomelos" <- 0,
  "Turdus_merula" <- 1.6,
  "Passer_domesticus" <- 3.6,
  "Sorex_araneus" <- 0.4,
  "Motacilla_cinerea" <- -3.6,
  "Motacilla_alba" <- -3.2,
  "Turdus_viscivorus" <- -3,
  "Fulica_atra" <- -2.6,
  "Cuculus_canorus" <- -2,
  "Apodemus_sylvaticus" <- 3.5,
  "Myodes_glareolus" <- 1.8,
  "Erithacus_rubecula" <- 1,
  "Sorex_minutus" <- 3.7,
  "Passer_domesticus" <- 2.2,
  
  "Carduelis_carduelis" <- -3.2,
  "Saxicola_rubetra" <- -3.4,
  "Pyrrhula_pyrrhula" <- -3.4,
  "Haematopus_ostralegus" <- -3.6,
  "Pluvialis_apricaria" <- -2.4,
  "Arvicola_amphibius" <- 2.3,
  "Periparus_ater" <- -2.3,
  "Prunella_modularis" <- -1.5,
  "Cyanystes_caeruleus" <- 3.45,
  
  "Cervus_elaphus" <- 1.6,
  "Cervus_nippon" <- 1.4,
  "Capreolus_capreolus" <- 1.8,
  "Troglodytes_troglodytes" <- 0.6,
  "Sciurus_carolinensis" <- -1.1,
  "Rana_temporaria" <- -0.3,
  "Aegithalos_caudatus" <- 1.2,
  "Loxia_curvirostra" <- -4.25,
  "Neomys_fodiens" <- 2.1,
  "Pipistrellus_pygmaeus" <- -2,
  "Lepus_timidus" <- 4,
  "Regulus_regulus" <- 3,
  "Acanthis_cabaret" <- -0.5,
  "Scolopax_rusticola" <- -1.3,
  "Pipistrellus_pipistrellus" <- -1.7,
  "Columba" <- -0.7,
  "Bufo_bufo" <- 0.3,
  "Anthus_trivialis" <- -2.7,
  "Fringilla_coelebs" <- -1.1,
  "Lophophanes_cristatus" <- -2.9,
  "Anatidae" <- 3.4,
  "Phylloscopus_trochilus" <- 3.2,
  "Tringa_nebularia" <- -4,
  "Actitis_hypoleucos" <- -1.5,
  "Certhia_familiaris" <- -4.4,
  "Sylvia_atricapilla" <- -0.85
  
)

fixed_x_positions[names(manual_positions)] <- manual_positions

# Set global boundaries
global_xlim <- c(min(fixed_x_positions)-0.2, max(fixed_x_positions)+0.2)
global_ylim <- c(1, 4)

## Initialize plot ---------------------------------------------------------
plot_filename <- "combined_trophic_web_aud.png"
png(plot_filename, width <- 2500, height <- 2500, res <- 600)

par(mar <- c(3,4,3,3), xaxs <- "i", yaxs <- "i")
plot.new()
plot.window(xlim <- global_xlim, ylim <- global_ylim)

## Plotting loop -----------------------------------------------------------
method_order <- c("Cam", "W", "T", "S", "Aud")

# Modify the plotting loop section like this:
for (method in method_order) {
  method_samples <- ncl_meta$sample[ncl_meta$type == method]
  detected_species <- colnames(ncl_pa)[colSums(ncl_pa[method_samples, ]) > 0]
  
  # Co-occurrence validation with explicit debugging
  method_links <- trophic_links %>%
    filter(resource %in% detected_species & 
             consumer %in% detected_species) %>%
    rowwise() %>%
    mutate(cooccurrence <- any(ncl_pa[method_samples, resource] == 1 & 
                                ncl_pa[method_samples, consumer] == 1)) %>%
    ungroup() %>%
    filter(cooccurrence) %>%
    filter(resource != consumer)
  
  # Add verification step
  cat("\nMethod:", method, "Links:", nrow(method_links), "\n")
  print(head(method_links))
  
  # Rest of your plotting code remains the same
  vertices_df <- data.frame(name <- detected_species) %>%
    left_join(link_pos, by <- c("name" <- "species")) %>%
    mutate(
      PreyAveragedTL <- as.numeric(PreyAveragedTL),
      x_position <- fixed_x_positions[match(name, names(fixed_x_positions))]
    )
  
  g <- graph_from_data_frame(
    d <- method_links[, c("resource", "consumer")],
    directed <- TRUE,
    vertices <- vertices_df
  )
  
  # Layout with fixed positions
  layout_matrix <- layout_with_fr(g, niter <- 500)
  layout_matrix[, 1] <- vertices_df$x_position
  layout_matrix[, 2] <- vertices_df$PreyAveragedTL
  
  cat("\nGraph edges for", method, ":")
  print(E(g))
  
  # Plot with NO LABELS
  plot.igraph(
    g,
    layout <- layout_matrix,
    vertex.size <- 0,
    vertex.color <- NA,
    vertex.frame.color <- NA,
    edge.arrow.size <- 0.2,
    edge.color <- adjustcolor(method_colors[method], alpha.f <- 0.4),
    edge.width <- 0.7,
    vertex.label <- NA,
    rescale <- FALSE,
    add <- TRUE,
    xlim <- global_xlim,
    ylim <- global_ylim
  )
  
  # # Add vertical labels manually - useful for initial understanding but removed from final plot
  # connected_nodes <- V(g)[degree(g, mode <- "all") > 0]
  # connected_layout <- layout_matrix[degree(g, mode <- "all") > 0, ]
  # 
  # text(
  #   x <- connected_layout[,1],
  #   y <- connected_layout[,2],
  #   labels <- connected_nodes$name,
  #   cex <- 0.4,
  #   srt <- 90,
  #   adj <- c(0, 0.5),
  #   xpd <- TRUE
  # )
}

## Final elements ----------------------------------------------------------
axis(2, at <- 1:4, las <- 1, cex.axis <- 0.5)
rect(
  xleft <- par("usr")[1],
  xright <- par("usr")[2],
  ybottom <- 1,
  ytop <- 4,
  border <- "black",
  lwd <- 2
)
mtext("Trophic Level", side <- 2, line <- 3, cex.lab <- 0.5)

legend(
  "bottom",
  legend <- names(method_colors),
  fill <- adjustcolor(method_colors, alpha.f <- 0.2),
  title <- "Methods",
  bty <- "n",
  horiz <- TRUE,
  cex <- 0.7,
  inset <- c(0.05, -0.15),
  xpd <- TRUE
)

dev.off()

#_____________________________________________
# Pooling for each method ----
library(dplyr)
library(igraph)
library(cheddar)
library(purrr)
library(combinat)
library(tidyr)
library(ggplot2)
library(stringr)

# Function to generate a co-occurrence network for one or more methods
generate_combined_network <- function(method_subset, ncl_pa, ncl_meta, trophic_links, link_pos) {
  # Get all sample names for the selected methods
  samples <- ncl_meta %>%
    filter(type %in% method_subset) %>%
    pull(sample)
  
  all_links <- list()
  all_species <- c()
  
  # Loop through each sample
  for (s in samples) {
    detected_species <- colnames(ncl_pa)[ncl_pa[s, ] == 1]
    
    # Subset trophic links for co-occurring species in that sample
    sample_links <- trophic_links %>%
      filter(resource %in% detected_species & consumer %in% detected_species)
    
    all_links[[s]] <- sample_links
    all_species <- union(all_species, detected_species)
  }
  
  # Combine all sample-specific links (de-duplicated)
  combined_links <- bind_rows(all_links) %>%
    distinct()  # remove duplicate rows across samples
  
  # Create vertices data frame for graph
  vertices_df <- data.frame(name <- all_species)
  
  # Construct igraph object
  g <- graph_from_data_frame(combined_links, directed <- TRUE, vertices <- vertices_df)
  
  return(list(graph <- g, species <- all_species, links <- combined_links))
}

# 1. Calculate metrics (without edge_density/fragmentation)
calculate_network_metrics <- function(g, detected_species, link_pos, method_name) {
  if (gorder(g) == 0 || gsize(g) == 0) {
    return(data.frame(
      type <- method_name,
      linkage_density <- NA,
      connectance <- NA,
      number_of_nodes <- NA,
      number_of_links <- NA,
      robustness <- NA,
      mean_trophic_level <- NA,
      modularity <- NA,
      avg_path_length <- NA
    ))
  }
  
  # for connectance calculation
  all_possible_links <- trophic_links %>% 
    distinct(resource, consumer) %>% 
    filter(resource %in% detected_species & consumer %in% detected_species)
  
  # Metrics
  linkage_density <- mean(degree(g))
  connectance <- if(nrow(all_possible_links) > 0) {
    gsize(g) / nrow(all_possible_links)
  } else {0}
  number_of_nodes <- gorder(g)
  number_of_links <- gsize(g)
  robustness <- calculate_robustness(g)
  
  trophic_levels <- link_pos %>%
    filter(species %in% detected_species) %>%
    pull(PreyAveragedTL)
  mean_trophic_level <- if (length(trophic_levels) > 0) mean(trophic_levels, na.rm <- TRUE) else NA
  
  clusters <- cluster_louvain(as_undirected(g))
  modularity <- modularity(clusters)
  avg_path_length <- mean_distance(g, directed <- FALSE)
  
  data.frame(
    type <- method_name,
    linkage_density <- linkage_density,
    connectance <- connectance,
    number_of_nodes <- number_of_nodes,
    number_of_links <- number_of_links,
    robustness <- robustness,
    mean_trophic_level <- mean_trophic_level,
    modularity <- modularity,
    avg_path_length <- avg_path_length
  )
}

# Generate network metrics
method_types <- unique(ncl_meta$type)
combo_metrics <- list()

for (i in 1:length(method_types)) {
  combos <- combn(method_types, i, simplify <- FALSE)
  for (combo in combos) {
    combo_name <- paste(combo, collapse <- "+")
    network_data <- generate_combined_network(combo, ncl_pa, ncl_meta, trophic_links, link_pos)
    
    metrics <- calculate_network_metrics(
      g <- network_data$graph,
      detected_species <- network_data$species,
      link_pos <- link_pos,
      method_name <- combo_name
    )
    
    combo_metrics[[combo_name]] <- metrics
  }
}

# Create filtered_metrics
all_combined_metrics <- bind_rows(combo_metrics)
filtered_metrics <- all_combined_metrics |> filter(!grepl("\\bP\\b", type))
filtered_metrics$perc_links <- (filtered_metrics$number_of_links / 183) * 100

# Generate metrics_long
metrics_long <- filtered_metrics |>
  pivot_longer(
    cols <- -c(type, number_of_nodes),
    names_to <- "metric",
    values_to <- "value"
  ) |> 
  mutate(num_methods <- str_count(type, "\\+") + 1)



# Plotting --------------
library(ggplot2)
library(cowplot)
library(stringr)
library(dplyr)

metrics <- c(
  "linkage_density", "perc_links", "connectance", "robustness", 
  "modularity", "mean_trophic_level"
)
axes_titles <- c(
  "Linkage Density", "% of Links", "Connectance", "Robustness (R50)", 
  "Modularity", "Mean Trophic Level"
)
title_mapping <- setNames(axes_titles, metrics)

highlight_combo <- "Aud+Cam+S+T+W"

# Color setup
method_colors <- c(
  "Aud" <- "#8de0da",
  "Cam" <- "#ee881c", 
  "S" <- "#f5c905",
  "T" <- "#2cd486",
  "W" <- "#435bf7",
  "Combined" <- "black"
)

method_labels <- c("All", "Aud", "Cam", "S", "T", "W")

# Plotting function
plot_metric_bar <- function(metric_name, data) {
  data_filtered <- data |>
    filter(metric == metric_name) |>
    filter(str_count(type, "\\+") == 0 | type == highlight_combo) |>
    mutate(
      type_renamed <- factor(
        ifelse(type == highlight_combo, "All", type),
        levels <- method_labels
      )
    )
  
  p <- ggplot(data_filtered, aes(x <- type_renamed, y <- value, fill <- type_renamed)) +
    geom_bar(stat <- "identity", alpha <- 0.6, width <- 0.8) +
    scale_fill_manual(values <- method_colors, guide <- "none") +
    labs(x <- NULL, y <- title_mapping[metric_name]) +
    theme_minimal() +
    theme(
      panel.background <- element_blank(),
      axis.line <- element_line(color <- "black"),
      axis.text.x <- element_text(angle <- 0, hjust <- 0.5, size <- 8),
      axis.text.y <- element_text(size <- 10),
      panel.grid <- element_blank()
    )
  
  # Adjust y-axis for mean trophic level
  if (metric_name == "mean_trophic_level") {
    max_y <- max(data_filtered$value, na.rm <- TRUE)
    p <- p + coord_cartesian(ylim <- c(2, max_y)) +
      scale_y_continuous(expand <- expansion(mult <- c(0, 0.05)))
  } else {
    p <- p + scale_y_continuous(expand <- expansion(mult <- c(0, 0.05)))
  }
  
  return(p)
}


# Generate and combine plots
metric_plots <- lapply(metrics, function(m) {
  plot_metric_bar(m, metrics_long)
})

final_plot <- plot_grid(
  plotlist <- metric_plots,
  ncol <- 2,
  labels <- paste0(LETTERS[4:(3 + length(metrics))], ")"),
  label_size <- 12,
  label_fontface <- "plain",
  label_x <- 0,
  label_y <- 1.05,
  hjust <- 0,
  align <- "v",
  axis <- "tblr"
) + theme(plot.margin <- unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

print(final_plot)

ggsave("final_metrics.png", final_plot)



# plotting all the combos

plot_metric_bar_all <- function(metric_name, data) {
  data_filtered <- data |>
    filter(metric == metric_name) |>
    mutate(
      num_methods <- str_count(type, "\\+") + 1,  # FIXED REGEX
      type_renamed <- ifelse(type == highlight_combo, "All", type)
    ) |>
    arrange(num_methods) |>  
    mutate(
      type_renamed <- factor(type_renamed, levels <- unique(type_renamed))
    )
  
  p <- ggplot(data_filtered, aes(x <- type_renamed, y <- value, fill <- type_renamed)) +
    geom_bar(stat <- "identity", alpha <- 0.7, width <- 0.8) +
    scale_fill_manual(values <- c(method_colors, "All" <- "black"), guide <- "none") +
    labs(x <- NULL, y <- title_mapping[metric_name]) +
    theme_minimal() +
    theme(
      panel.background <- element_blank(),
      axis.line <- element_line(color <- "black"),
      axis.text.x <- element_text(angle <- 45, hjust <- 1, size <- 8),
      axis.text.y <- element_text(size <- 10),
      panel.grid <- element_blank()
    ) +
    scale_y_continuous(expand <- expansion(mult <- c(0, 0.05)))
  
  if (metric_name == "mean_trophic_level") {
    max_y <- max(data_filtered$value, na.rm <- TRUE)
    p <- p + coord_cartesian(ylim <- c(2, max_y))
  }
  
  return(p)
}


metric_plots <- lapply(metrics, function(m) {
  plot_metric_bar_all(m, metrics_long)
})

final_plot <- plot_grid(
  plotlist <- metric_plots,
  ncol <- 2,
  labels <- paste0(LETTERS[4:(3 + length(metrics))], ")"),
  label_size <- 12,
  label_fontface <- "plain",
  label_x <- 0,
  label_y <- 1.05,
  hjust <- 0,
  align <- "v",
  axis <- "tblr"
) + theme(plot.margin <- unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

print(final_plot)

# With a bar for the entire potential network ----------------------------------
# 1. Calculate NCL Metrics (unchanged)
calculate_connectance <- function(num_links, possible_links) {
  ifelse(possible_links > 0, num_links / possible_links, NA)
}

NCL.num_links <- NumberOfTrophicLinks(NCL.community)
NCL.linkage_density <- NCL.num_links / NumberOfNodes(NCL.community)
possible_links <- NumberOfNodes(NCL.community) * (NumberOfNodes(NCL.community) - 1)
NCL.connectance <- calculate_connectance(NCL.num_links, possible_links)
graph_for_robustness <- graph_from_data_frame(trophic_links, directed <- TRUE, vertices <- data.frame(name <- nodes$node))
NCL.robustness <- calculate_robustness(graph_for_robustness)
NCL.mean_trophic_level <- mean(link_pos$PreyAveragedTL)

# 2. Create 'All potential' baseline
NCL_all_potential <- data.frame(
  type <- "All potential",
  linkage_density <- NCL.linkage_density,
  connectance <- NCL.connectance,
  robustness <- NCL.robustness,
  mean_trophic_level <- NCL.mean_trophic_level
)

# 3. Combine with existing data
metrics_long <- bind_rows(
  NCL_all_potential |> pivot_longer(-type, names_to <- "metric", values_to <- "value"),
  metrics_long  # Your original data containing other methods
)

# 4. Updated plotting parameters
metrics <- c("linkage_density", "connectance", "robustness", "mean_trophic_level")
axes_titles <- c("Linkage Density", "Connectance", "Robustness (R50)", "Mean Trophic Level")
title_mapping <- setNames(axes_titles, metrics)
highlight_combo <- "Cam+Aud+S+T+W"

method_colors <- c(
  "All potential" <- "black",
  "Aud" <- "#8de0da",
  "Cam" <- "#ee881c", 
  "S" <- "#f5c905",
  "T" <- "#2cd486",
  "W" <- "#435bf7",
  "Combined" <- "darkgrey"
)

method_labels <- c("All potential", "Combined", "Aud", "Cam", "S", "T", "W")

# 5. Modified plotting function
plot_metric_bar <- function(metric_name, data) {
  data_filtered <- data |>
    filter(metric == metric_name) |>
    filter(type == "All potential" | 
             str_count(type, "\\+") == 0 | 
             type == highlight_combo) |>
    mutate(
      type_renamed <- factor(
        case_when(
          type == "All potential" ~ "All potential",
          type == highlight_combo ~ "Combined",
          TRUE ~ type
        ),
        levels <- method_labels
      )
    )
  
  p <- ggplot(data_filtered, aes(x <- type_renamed, y <- value, fill <- type_renamed)) +
    geom_bar(stat <- "identity", alpha <- 0.7, width <- 0.8) +
    scale_fill_manual(values <- method_colors, guide <- "none") +
    labs(x <- NULL, y <- title_mapping[metric_name]) +
    theme_minimal() +
    theme(
      panel.background <- element_blank(),
      axis.line <- element_line(color <- "black"),
      axis.text.x <- element_text(angle <- 0, hjust <- 0.5, size <- 8),
      axis.text.y <- element_text(size <- 10),
      panel.grid <- element_blank()
    )
  
  if (metric_name == "mean_trophic_level") {
    max_y <- max(data_filtered$value, na.rm <- TRUE)
    p <- p + coord_cartesian(ylim <- c(2, max_y)) +
      scale_y_continuous(expand <- expansion(mult <- c(0, 0.05)))
  } else {
    p <- p + scale_y_continuous(expand <- expansion(mult <- c(0, 0.05)))
  }
  
  return(p)
}

# 6. Generate final plot (unchanged)
metric_plots <- lapply(metrics, function(m) {
  if (!m %in% metrics_long$metric) warning(paste(m, "not found in data!"))
  plot_metric_bar(m, metrics_long)
})

final_plot <- plot_grid(
  plotlist <- metric_plots,
  ncol <- 2,
  labels <- paste0(LETTERS[4:(3 + length(metrics))], ")"),
  label_size <- 12,
  label_fontface <- "plain",
  label_x <- 0,
  label_y <- 1.05,
  hjust <- 0,
  align <- "v",
  axis <- "tblr"
) + theme(plot.margin <- unit(c(0.5, 0.5, 0.5, 0.5), "cm"))

print(final_plot)
