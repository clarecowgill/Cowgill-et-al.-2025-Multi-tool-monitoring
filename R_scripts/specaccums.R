
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(iNEXT)
library(lubridate)

#______________________________________________________________________________
# Load data ----

# transpose filtered read count matrix
ncl_sac<-as.data.frame(read.csv('ncl_all_pa.csv', row.names <- 1))
nclvegmet <- read.csv("all_meta.csv")


# merge by sample to bring in habitat and type data
ncl_sac$sample <- row.names(ncl_sac)
ncl_sac <- merge(ncl_sac, nclvegmet, by <- "sample")

# ______________________________________________________________________________
# iNEXT accumulation curves ----
## SAC split by sample type ----

# Function to create presence/absence matrix for a given type
pa_matrix_type <- function(df, type_value) {
  # Filter by the specified type
  filtered_df <- df %>% filter(type == type_value)

  rownames(filtered_df)<-filtered_df$sample
  
  numeric_cols<-filtered_df %>% select_if(is.numeric)
  
  pa_matrix<-numeric_cols
  pa_matrix[pa_matrix > 0]<-1
  pa_matrix[pa_matrix <= 0]<-0
  
  rownames(pa_matrix)<-rownames(filtered_df)
  
  return(t(pa_matrix))
}

# Create matrices
matrix_T<-pa_matrix_type(ncl_sac, "T")
matrix_W<-pa_matrix_type(ncl_sac, "W")
matrix_S<-pa_matrix_type(ncl_sac, "S")
matrix_P<-pa_matrix_type(ncl_sac, "P")

# Combine matrices into a list
type_matrix_list<-list(
  'Tree roller'<-matrix_T,
  'Water'<-matrix_W,
  'Soil'<-matrix_S,
  'Scat'<-matrix_P
)

# add all 

sac_dna<-ncl_sac %>% 
  filter(type %in% c('T', 'S', 'P', 'W'))

# Function to create a combined presence/absence matrix for all types
pa_matrix_all<-function(df) {
  rownames(df)<-df$sample
  
  # Select only numeric columns (species presence/absence data)
  numeric_cols<-df %>% select_if(is.numeric)
  
  # Convert to presence/absence
  pa_matrix<-numeric_cols
  pa_matrix[pa_matrix > 0]<-1
  pa_matrix[pa_matrix <= 0]<-0
  
  rownames(pa_matrix)<-df$sample
  
  return(t(pa_matrix))  # Transpose to match the other matrices
}

# Create matrix for all types combined
matrix_All<-pa_matrix_all(sac_dna)

# Combine matrices into a list, including "All"
type_matrix_list<-list(
  'Tree roller'<-matrix_T,
  'Water'<-matrix_W,
  'Soil'<-matrix_S,
  'Scat'<-matrix_P,
  'All'<-matrix_All  # Adding the combined matrix
)

# Apply iNEXT to each matrix in the list
#Chao1, 999 permutations, bootstrap 50 - extrapolation
unpooled_250_output<-iNEXT(type_matrix_list, datatype<-"incidence_raw", 
                            conf<-0.95, se<-TRUE, endpoint<-250)

# Update the plot with custom colors and labels, including "All"
unpooled_250<-ggiNEXT(unpooled_250_output) +
  scale_colour_manual(values=c("#555555","maroon", "#edd134", "#349945", "#435bf7")) +  # Adding color for "All"
  scale_fill_manual(values=c("#555555","maroon", "#edd134", "#349945", "#435bf7")) +
  labs(x<-'Sample number') +
  labs(y<-'Taxonomic richness') +
  ylim(0, 85) + # Set the y-axis limit to go from 0 to 100 to match other plot
  theme_minimal() +
  theme(
    axis.line<-element_line(color<-"black"),
    panel.grid.major<-element_blank(),
    panel.grid.minor<-element_blank(),
    panel.background<-element_blank(),
    plot.background<-element_blank(),
    legend.position<-"bottom",
    axis.title.x<-element_text(size<-15),
    axis.title.y<-element_text(size<-15),
    axis.text.x<-element_text(size<-12),
    axis.text.y<-element_text(size<-12)
  )

ggsave("sac_extrap_unpooled_250_all.png", plot<-unpooled_250, dpi<-800, width<-4.5, height<-5.2)

# ______________________________________________________________________________
## Method comparison with camera and acoustic data ----

### Camera and Audiomoth deployment days total ----

ncl_cams<-read.csv('ncl_cams.csv')
ncl_aud<-read.csv('ncl_acoust.csv')

# function to create presence/absence matrix
method_pa_matrix<-function(df, method_value, group_by_cols, expand_cols) {
  # Define fixed date ranges for each Type
  date_ranges<-list(
    Acoustic<-seq(as.Date("2023-05-01"), as.Date("2023-06-15"), by<-"day"),
    Camera<-seq(as.Date("2023-05-01"), as.Date("2023-06-30"), by<-"day")
  )
  
  # Convert the date column to date object
  df$date<-dmy(df$date)
  
  # Get the relevant date range based on the method
  fixed_dates<-date_ranges[[method_value]]
  
  # Get lists of unique species and sites
  all_species<-unique(df$species)
  all_sites<-unique(df$site)
  
  # Create a complete grid of combinations for dates, sites, and species
  complete_grid<-expand.grid(date<-fixed_dates, site<-all_sites, 
                              species<-all_species)
  
  # Summarise to presence/absence
  presence_absence<-df %>%
    group_by(across(all_of(group_by_cols))) %>%
    summarise(present<-ifelse(n() > 0, 1, 0), .groups<-'drop') %>%
    right_join(complete_grid, by<-c("date", "site", "species")) %>%
    mutate(present<-ifelse(is.na(present), 0, present)) %>%
    pivot_wider(names_from<-all_of(expand_cols), values_from<-present, 
                values_fill<-list(present<-0))
  
  # Convert to numeric matrix
  pa_matrix<-as.matrix(presence_absence %>% dplyr::select(., -species))
  rownames(pa_matrix)<-presence_absence$species
  
  return(pa_matrix)
}

# Create matrices for camera and acoustic data
matrix_cam_days<-method_pa_matrix(ncl_cams, "Camera", 
                                   c("date", "site", "species"), 
                                   c("date", "site"))
matrix_aud_days<-method_pa_matrix(ncl_aud, "Acoustic", 
                                   c("date", "site", "species"), 
                                   c("date", "site"))

# Harmonize species (rows)
all_species<-union(rownames(matrix_cam_days), rownames(matrix_aud_days))

pad_rows<-function(mat, all_species) {
  missing_species<-setdiff(all_species, rownames(mat))
  if (length(missing_species) > 0) {
    missing_mat<-matrix(0, nrow<-length(missing_species), ncol<-ncol(mat))
    rownames(missing_mat)<-missing_species
    colnames(missing_mat)<-colnames(mat)
    mat<-rbind(mat, missing_mat)
  }
  return(mat[all_species, , drop<-FALSE])
}

matrix_cam_days_full<-pad_rows(matrix_cam_days, all_species)
matrix_aud_days_full<-pad_rows(matrix_aud_days, all_species)

# Harmonize date-site combinations (columns)
all_cols<-union(colnames(matrix_cam_days_full), colnames(matrix_aud_days_full))

pad_cols<-function(mat, all_cols) {
  missing_cols<-setdiff(all_cols, colnames(mat))
  if (length(missing_cols) > 0) {
    missing_mat<-matrix(0, nrow<-nrow(mat), ncol<-length(missing_cols))
    rownames(missing_mat)<-rownames(mat)
    colnames(missing_mat)<-missing_cols
    mat<-cbind(mat, missing_mat)
  }
  return(mat[, all_cols, drop<-FALSE])
}

matrix_cam_days_full<-pad_cols(matrix_cam_days_full, all_cols)
matrix_aud_days_full<-pad_cols(matrix_aud_days_full, all_cols)

# Combine matrices
matrix_combined_days<-matrix_cam_days_full + matrix_aud_days_full
matrix_combined_days[matrix_combined_days > 1]<-1  # Make binary

# List for iNEXT
cam_aud_list_days<-list(
  'Acoustic'<-matrix_aud_days,
  'Camera'<-matrix_cam_days,
  'Both'<-matrix_combined_days
)

# Run iNEXT
days_output<-iNEXT(cam_aud_list_days, datatype<-"incidence_raw", 
                    conf<-0.95, se<-TRUE, endpoint<-1500)

# Plot with three lines
days_output_plot<-ggiNEXT(days_output) +
  scale_colour_manual(
    values<-c(
      'Acoustic'<-"#8de0da", 'Camera'<-"#ee881c", 'Both'<-"#555555")) +
  scale_fill_manual(
    values<-c('Acoustic'<-"#8de0da",'Camera'<-"#ee881c",'Both'<-"#555555")) +
  labs(x<-'Number of deployment days',
       y<-'Taxonomic richness') +
  ylim(0, 85) +
  theme_minimal() +
  theme(
    axis.line<-element_line(color<-"black"),
    panel.grid.major<-element_blank(),
    panel.grid.minor<-element_blank(),
    panel.background<-element_blank(),
    plot.background<-element_blank(),
    legend.position<-"bottom",
    axis.title.x<-element_text(size<-15),
    axis.title.y<-element_text(size<-15),
    axis.text.x<-element_text(size<-12),
    axis.text.y<-element_text(size<-12)
  )

# Show plot
days_output_plot

# Save plot
ggsave("cam_aud_combined_deploymentdays_sac.png", plot<-days_output_plot, 
       dpi<-300, width<-4.5, height<-5.2)

