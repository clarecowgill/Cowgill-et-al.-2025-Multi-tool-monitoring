

library(vegan)
library(ggplot2)
library(viridis)
library(ggpubr)
library(dplyr)
library(tidyr)
library(cowplot)
library(dunn.test)

#______________________________________________________________________________
# Load data ----

methods_week_pa <- read.csv('Data/all_week_pa.csv', row.names <- 1)
week_ca_vegmet <- read.csv('Data/all_week_meta.csv', row.names <- 1)

# ______________________________________________________________________________

# ALPHA DIVERSITY ----
method_colors <- c(
    "Aud" <- "#8de0da",
    "Cam" <- "#ee881c", 
    "P" <- 'maroon',
    "S" <- "#edd134",
    "T" <- "#349945",
    "W" <- "#435bf7"
)


# Diversity plot function
plot_diversity <- function(data, index, index_label) {
  p <- ggplot(data, aes_string(x <- "type", y <- index, fill <- "type")) + 
    geom_boxplot(alpha <- 0.7, color <- "black") +
    geom_jitter(colour <- "black", alpha <- 0.4, position <- position_jitter(width <- 0.1), size <- 0.8) +
    scale_x_discrete(labels <- c("P" <- "Scat", "S" <- "Soil", "T" <- "Tree", "W" <- "Water")) +
    scale_fill_manual(values <- method_colors) +
    ylab(index_label) + 
    xlab(NULL) +
    theme_minimal(base_size <- 15) +
    theme(
      axis.line <- element_line(color <- "black"),
      panel.grid.major <- element_blank(),
      panel.grid.minor <- element_blank(),
      panel.background <- element_blank(),
      plot.background <- element_blank(),
      legend.position <- "none"
    )
  return(p)
}

# Filter week_abundance
filtered_samples <- week_ca_vegmet$sample
methods_week <- methods_week_pa[rownames(methods_week_pa) %in% filtered_samples, ]

# Calculate diversity metrics
methods_week$richness <- specnumber(methods_week)

# Add sample column and merge
methods_week$sample <- rownames(methods_week)
merged_data <- merge(methods_week, week_ca_vegmet, by <- "sample")

richness_plot <- plot_diversity(merged_data, "richness", "Taxonomic richness")
richness_plot

ggsave("richness_boxplot.png", plot <- richness_plot, width <- 5, height <- 5,  dpi <- 800, bg <- 'white')

dunn.test(merged_data$richness, merged_data$type, method <- "bh")

median_richness <- merged_data %>%
  group_by(type) %>%
  summarize(median_richness <- median(richness))

# ______________________________________________________________________________
# BETA DIVERSITY ----

# By method type, expanding acoustics and cams by day/week ----

## Jaccard only for presence/ absence

all_ncl_jac_week <- vegdist(methods_week_pa, method <- "jaccard")

# Reorder meta data to match
week_ca_vegmet <- week_ca_vegmet[order(week_ca_vegmet$sample), ]

# Compute homogeneity of group dispersions (variances)
all_jac_var_week <- betadisper(all_ncl_jac_week, week_ca_vegmet$type)

# Check homogeneity of multivariate dispersions. Groups being tested 
# should have the sample multivariate spread to conform to assumptions of PERMANOVA.
all_mod_jac_week <- with(week_ca_vegmet, all_jac_var_week)

# Compute mean distance to centroid and variance per group
tapply(all_jac_var_week$distances, week_ca_vegmet$type, mean)
tapply(all_jac_var_week$distances, week_ca_vegmet$type, var)

# ordination plot of distance to centroids
plot(all_jac_var_week)

# boxplot of distance to centroids
boxplot(all_jac_var_week, xlab<-"Type", xaxt<-"n", bty<-"n", main<-"Boxplot of all_jac_var_week")
axis(side<-1, at<-1:6, labels<-c("Aud", "Cam", "Tree", "Scat", "Soil", "Water"))

#check whether variance is different between locations using standard parametric anova or permutation tests
anova(all_jac_var_week)
permutest(all_jac_var_week)

# Analyse pairwise differences between locations using parametric Tukey's HSD test
TukeyHSD(all_jac_var_week)

# ordination of beta diversity (obd) with metaMDS
all_obd_jac_week <- metaMDS(all_ncl_jac_week,
                      dist <- "jaccard",
                      k<-3,
                      maxit<-999,
                      trymax<-1000,
                      noshare<-TRUE,
                      wascores<-TRUE)

# Assess goodness of ordination fit (stress plot)
all_obd_jac_week$stress
stressplot(all_obd_jac_week)

## Plot site scores as text
ordiplot(all_obd_jac_week, display <- "sites", type <- "text", cex<-0.5)


## Build dataframes with NMDS coordinates and metadata
all_ncl_jac_NMDS1 <- all_obd_jac_week$points[,1]
all_ncl_jac_NMDS2 <- all_obd_jac_week$points[,2]
all_ncl_jac_NMDS <- data.frame(NMDS1<-all_ncl_jac_NMDS1, 
                              NMDS2<-all_ncl_jac_NMDS2,
                              type <- week_ca_vegmet$type)


# Create a new column 'type_grouped' where 'P', 'T', 'S', 'W' are combined into one group
all_ncl_jac_NMDS$type_grouped <- ifelse(all_ncl_jac_NMDS$type %in% c("P", "T", "S", "W"), "DNA", all_ncl_jac_NMDS$type)

# Ensure Grouped rows are first in the dataset
all_ncl_jac_NMDS$type_grouped <- factor(all_ncl_jac_NMDS$type_grouped, 
                                        levels <- c("DNA", "Aud", "Cam"))
# Reorder factor levels so "Cam" comes last
all_ncl_jac_NMDS$type <- factor(all_ncl_jac_NMDS$type,
                                levels <- c("Aud", "P", "S", "T", "W", "Cam")
)

# switch type and type grouped to change grouping 
all_ncljacnmds_plot <- ggplot(all_ncl_jac_NMDS, aes(x <- NMDS1, y <- NMDS2, colour <- type, fill <- type)) + 
  geom_point(cex <- 4, alpha <- 0.9) + 
  stat_ellipse(alpha<-0.3, geom<-"polygon", aes(fill<-type), colour<-NA) +
  labs(x <- "NMDS1", y <- "NMDS2", subtitle <- paste("stress <-", round(all_obd_jac_week$stress, 3))) +  
  scale_colour_manual(name <- "type",
                      values <- c("#8de0da","maroon","#edd134", '#349945', '#435bf7', '#ee881c')) +
  scale_fill_manual(name <- "type",
                    values <- c("#8de0da","maroon","#edd134", '#349945', '#435bf7', '#ee881c')) +
  coord_cartesian(xlim <- c(-3.3, 2.7), ylim <- c(-3.8, 3)) + 
  theme(panel.background <- element_rect(fill <- 'white'),
        axis.line.x <- element_line(colour <- 'black', linewidth <- 0.5, linetype<-'solid'),
        axis.line.y <- element_line(colour <- 'black', linewidth <- 0.5, linetype<-'solid'),
        axis.title.x <- element_text(margin <- unit(c(8, 0, 0, 0), "mm")),
        axis.title.y <- element_text(margin <- unit(c(0, 5, 0, 0), "mm")),
        axis.text.x <- element_text(colour<-"black"),
        axis.text.y <- element_text(colour<-"black"),
        plot.title <- element_text(face<-"bold", hjust<-0, colour<-"black"),
        plot.subtitle <- element_text(face <- "plain", hjust <- 0, color <- "black",
                                     margin <- unit(c(2, 0, 0, 0), "mm"), size <- 18),
        text <- element_text(size<-14),
        legend.position <- "bottom",
        legend.direction<-"horizontal",
        legend.box<-"vertical",
        legend.key<-element_blank())

all_ncljacnmds_plot

ggsave(method_jac_nmds_week_all.png', plot<- all_ncljacnmds_plot, width <- 8, height <- 8, dpi <- 300)

# Create a new column 'type_grouped' where 'P', 'T', 'S', 'W' are combined into one group for meta data
week_ca_vegmet$type_grouped <- ifelse(week_ca_vegmet$type %in% c("P", "T", "S", "W"), "DNA", week_ca_vegmet$type)

# Run permanova with adonis to check for difference in spatial turnover of communities
week_jac_adonis <- adonis2(all_ncl_jac_week ~ type, week_ca_vegmet)

week_jac_adonis


# check for other differences
dispersion_test <- betadisper(all_ncl_jac_week, week_ca_vegmet$type)
permutest(dispersion_test, permutations <- 999)

TukeyHSD(dispersion_test)

anosim_result <- anosim(all_ncl_jac_week, week_ca_vegmet$type)
summary(anosim_result)


# 3D plot to explore the data better (since k<-3)
library(plotly)
all_ncl_jac_NMDS3 <- all_obd_jac_week$points[,3]  # Assuming the third dimension is available

# Create a data frame with the three dimensions and metadata
all_ncl_jac_NMDS_3D <- data.frame(NMDS1 <- all_ncl_jac_NMDS1, 
                                  NMDS2 <- all_ncl_jac_NMDS2, 
                                  NMDS3 <- all_ncl_jac_NMDS3,
                                  type <- week_ca_vegmet$type)

# 3D scatter plot
plot_ly(all_ncl_jac_NMDS_3D, 
        x <- ~NMDS1, 
        y <- ~NMDS2, 
        z <- ~NMDS3, 
        color <- ~type, 
        colors <- c("#05925b", "#ee881c", "maroon", "#edd134", 'darkgreen', 'blue'),
        type <- 'scatter3d', 
        mode <- 'markers', 
        marker <- list(size <- 5, opacity <- 0.8)) %>%
  layout(scene <- list(xaxis <- list(title <- 'NMDS1'),
                      yaxis <- list(title <- 'NMDS2'),
                      zaxis <- list(title <- 'NMDS3')),
         title <- paste("3D NMDS Plot (Stress <-", round(all_obd_jac_week$stress, 3), ")"))


# Load the package
library(pairwiseAdonis)

# Run pairwise PERMANOVA
pairwise_results <- pairwise.adonis(all_ncl_jac_week, factors <- week_ca_vegmet$type)

# View the results
pairwise_results

    