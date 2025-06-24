
library(lubridate)
library(dplyr)
library(data.table)
library(tidyr)
library(stringr)

# ______________________________________________________________________________

# handling blank contamination with LOC and LOQ ----

ncl <- read.csv('ncl_matrix_raw.csv')
blank_association <- read.csv('blank_association_all.csv')
type_dat <- read.csv('edna_meta.csv')

ncl_l <- as.data.frame(as.table(as.matrix(ncl)))
colnames(ncl_l) <- c("species", "sample", "reads")

# Merge blank_association with type data to associate blanks with sample types
blank_type_association <- merge(blank_association, type_dat, by <- "sample")

# Subset to only get the blank reads
blank_reads <- ncl_l[ncl_l$sample %in% blank_type_association$blank, ]
blank_reads$reads <- as.numeric(trimws(blank_reads$reads))

# Merge blank reads with the blank associations to link the reads with the sample types
merged_blank_data <- merge(blank_reads, blank_type_association, 
                          by.x <- "sample", by.y <- "blank", all.x <- TRUE)

merged_blank_data$reads <- as.numeric(trimws(merged_blank_data$reads))


# Calculate mean and standard deviation of reads for each species/sample_type combination
blank_stats <- aggregate(reads ~ species + type, data <- merged_blank_data, 
                        FUN <- function(x) c(mean <- mean(x), sd <- sd(x)))

# Convert to data frame and rename columns
blank_stats <- do.call(data.frame, blank_stats)
colnames(blank_stats) <- c("species", "type", "mean_blank_reads", "sd_blank_reads")

# Calculate LOD and LOQ
blank_stats$LOD <- blank_stats$mean_blank_reads + 3 * blank_stats$sd_blank_reads
blank_stats$LOQ <- blank_stats$mean_blank_reads + 10 * blank_stats$sd_blank_reads

# Merge blank_stats with the original ncl data based on species and sample type
non_blank_data <- merge(ncl_l, type_dat, by <- "sample")
lod_loq_data <- merge(non_blank_data, blank_stats, by <- c("species", "type"))
lod_loq_data$reads <- as.numeric(trimws(lod_loq_data$reads))

# filter data to change reads less than the LOD to 0
below_lod <- lod_loq_data$reads < lod_loq_data$LOD
lod_loq_data$reads[below_lod] <- 0

# turn filtered data back into a matrix
filtered_data_matrix <- lod_loq_data[, c("species", "sample", "reads")]

# Convert to wide format using reshape
wide_data <- reshape(filtered_data_matrix, 
                     idvar <- "species", 
                     timevar <- "sample", 
                     direction <- "wide")
colnames(wide_data) <- gsub("reads.", "", colnames(wide_data), fixed <- TRUE)
rownames(wide_data) <- wide_data$species
wide_data$species <- NULL

# add taxonomy
ncl$species <- rownames(ncl)
taxonomy <- unique(ncl[, c("species", "taxonomy")])
wide_data$species <- rownames(wide_data)
ncl_fil <- merge(wide_data, taxonomy, by <- "species", all.x <- TRUE)
rownames(ncl_fil) <- ncl_fil$species
ncl_fil$species <- NULL

# to use the species-specific filtered data
ncl <- ncl_fil
#______________________________________________________________________________

# filter data to remove blanks, fish and domestic species ----

# grab taxonomy:
ncl_tax <- ncl$taxonomy

# drop taxonomy:
ncl <- ncl[!names(ncl) %in% 'taxonomy']

# exclude irrelevant species:
excluded_taxa <- c('Actinopteri', 'Hyperoartia', 
                  'Homo_sapiens', 'Sus_scrofa', 
                  'Canis_lupus', 'Ovis_aries', 'Bos_taurus', 'Bovidae')
excluded <- paste(excluded_taxa, collapse <- '|')
ncl <- ncl[grep(excluded, ncl_tax, invert <- T),]

# drop controls:
blanks <- c('NEG', 'POS', 'EB', "BK", "FB")
blanks <- paste(blanks, collapse <- '|')
ncl <- ncl[grep(blanks, names(ncl), invert <- T)]

# cleanup option 1
# 0.1% threshold of total reads in a sample:
ncl[t(t(ncl)/colSums(ncl)) < 0.001] <- 0

write.csv(ncl, 'ncl_filtered.csv')


# _____________________________________________________________________________
# Filter ncl data to remove overlap across taxonomic levels

ncl_f <- ncl

combine_species <- function(data, species_to_combine, new_species_name) {
  # Subset the rows for the species to combine
  species_subset <- data[species_to_combine, , drop <- FALSE]
  
  # Ensure that only numeric values are summed
  # species_subset_numeric <- as.data.frame(lapply(species_subset, as.numeric))
  
  # Sum those rows across the samples (with NA values ignored)
  summed_data <- colSums(species_subset, na.rm <- TRUE)
  
  # Remove the original species rows
  data <- data[!rownames(data) %in% species_to_combine, ]
  
  # Add a new row in the dataframe with the summed values
  data[new_species_name, ] <- summed_data
  
  return(data)
}

ncl_f <- combine_species(ncl_f, c("Anas_carolinensis", "Anas_crecca", "Anatidae", "Aythya"), "Anatidae")
ncl_f <- combine_species(ncl_f, c("Cervidae", "Cervus_nippon"), "Cervus_nippon")
ncl_f <- combine_species(ncl_f, c("Columba", "Columba_livia", "Columba_oenas"), "Columba")
ncl_f <- combine_species(ncl_f, c("Phasianus_colchicus", "Phasianidae"), "Phasianidae")
ncl_f <- combine_species(ncl_f, c("Fringilla", "Fringillidae"), "Fringilla_coelebs")

rownames(ncl_f)[rownames(ncl_f) == "Regulus"] <- "Regulus_regulus"
rownames(ncl_f)[rownames(ncl_f) == "Passer"] <- "Passer_domesticus"
rownames(ncl_f)[rownames(ncl_f) == "Phylloscopus"] <- "Phylloscopus_trochilus"
rownames(ncl_f)[rownames(ncl_f) == "Turdus"] <- "Turdus_merula"

ncl_f <- ncl_f[rownames(ncl_f) != "Passeriformes", ]

write.csv(ncl_f, 'ncl_f.csv')

# _____________________________________________________________________________
library(reshape2)

# prep data for beta diversity analysis ----
nclvegmet <- read.csv('edna_meta.csv')
# presence/absence:
ncl_pa <- ncl_f
ncl_pa[ncl_pa > 0] <- 1
ncl_pa <- t(ncl_pa)

# species composition (proportion reads)
ncl_numeric <- ncl[, sapply(ncl, is.numeric)]
sapply(ncl_numeric, class)
ncl_prc <- t(ncl_numeric) / colSums(ncl_numeric)

# Drop samples with nowt in em
empty_samples <- ncl_pa[rowSums(ncl_pa) == 0, ]
ncl_pa <- ncl_pa[!rowSums(ncl_pa) == 0,]
ncl_prc <- ncl_prc[!is.na(rowSums(ncl_prc)) & rowSums(ncl_prc) != 0, ]

kept_samples <- rownames(ncl_pa)
nclvegmet_fil <- nclvegmet[nclvegmet$sample %in% kept_samples, ]

ncl_prc_dt <- as.data.table(ncl_prc, keep.rownames <- "sample")  # rownames become 'sample'
ncl_prc_long <- melt(ncl_prc_dt, id.vars <- "sample", variable.name <- "species", value.name <- "proportion")

write.csv(ncl_pa, 'ncl_pa.csv')
write.csv(ncl_prc, 'ncl_prc.csv')
write.csv(nclvegmet_fil, 'nclvegmet_fil.csv')


# ______________________________________________________________________________

# make a long version for bubble plots ----

# Convert row names to a column named 'species'
ncl_f$species <- rownames(ncl_f)

# use reshape to convert from wide to long format
ncl_long <- reshape(ncl_f,
                    varying <- list(names(ncl_f)[!names(ncl_f) %in% c("species")]),
                    v.names <- "reads", 
                    timevar <- "sample", 
                    times <- names(ncl_f)[!names(ncl_f) %in% c("species")],
                    idvar <- "species", 
                    direction <- "long")
rownames(ncl_long) <- NULL

write.csv(ncl_long,'ncl_long.csv')


# _____________________________________________________________________________

### data for beta diversity with cameras and acoustics ----

# just presence absence to begin - jaccard analysis

ncl_cams <- read.csv('ncl_cams.csv')
ncl_aud <- read.csv('ncl_acoust.csv')

ncl_aud$species <- gsub(" ", "_", ncl_aud$species)
ncl_cams$species <- gsub(" ", "_", ncl_cams$species)


ncl_camaud <- rbind(ncl_cams, ncl_aud)
unique(ncl_camaud$species)

# presence/absence:
ncl_ca_pa <- with(ncl_camaud, table(site, species))
ncl_ca_pa[ncl_ca_pa > 0] <- 1

### combining into a methods matrix ----

ncl_long <- ncl_long[ncl_long$reads != 0, ]
ncl_long <- ncl_long[, -3]

ncl_camaud <- ncl_camaud[, -3]
ncl_camaud <- ncl_camaud[, -4]
ncl_camaud_subset <- ncl_camaud[, -2]
names(ncl_camaud_subset)[names(ncl_camaud_subset) == "site"] <- "sample"

ncl_all <- rbind(ncl_camaud_subset, ncl_long)
unique(ncl_all$species)

# presence/absence:
ncl_all_pa <- with(ncl_all, table(sample, species))
ncl_all_pa[ncl_all_pa > 0] <- 1

write.csv(ncl_all_pa, 'ncl_all_pa_unfiltered.csv')

# filter dataframe to remove overlap and repeats of family or genus level data

ncl_all_pa <- as.data.frame.matrix(ncl_all_pa)

merge_columns <- function(data, col1, col2, new_col) {
  
  data[[new_col]] <- as.integer(data[[col1]] == 1 | data[[col2]] == 1)
  # Remove the old columns
  data[[col1]] <- NULL
  data[[col2]] <- NULL
  return(data)
  
}

colnames(ncl_all_pa)

ncl_all_pa$Parus_major <- ifelse(ncl_all_pa$Parus_major == 1 | ncl_all_pa$Paridae == 1, 1, 0)
ncl_all_pa[["Paridae"]] <- NULL

ncl_all_pa$Columba <- ifelse(ncl_all_pa$Columba == 1 | ncl_all_pa$Columba_oenas == 1, 1, 0)
ncl_all_pa[["Columba_oenas"]] <- NULL

ncl_all_pa$Columba <- ifelse(ncl_all_pa$Columba == 1 | ncl_all_pa$Columba_palumbus == 1, 1, 0)
ncl_all_pa[["Columba_palumbus"]] <- NULL

ncl_all_pa$Anatidae <- ifelse(ncl_all_pa$Anatidae == 1 | ncl_all_pa$Anas_crecca == 1, 1, 0)
ncl_all_pa[["Anas_crecca"]] <- NULL

ncl_all_pa$Phasianidae <- ifelse(ncl_all_pa$Phasianidae == 1 | ncl_all_pa$Phasianus_colchicus == 1, 1, 0)
ncl_all_pa[["Phasianus_colchicus"]] <- NULL

ncl_all_pa[["Passeriformes"]] <- NULL # too vague for current analysis

write.csv(ncl_all_pa, 'ncl_all_pa.csv')


# _____________________________________________________________________________
## by-week data for cams and acoustics -----
library(lubridate)

ca_vegmet <- read.csv("all_meta.csv", row.names <- 1)
ca_vegmet$site <- row.names(ca_vegmet)

# Ensure the date column in cams is in Date format and add week identifier
ncl_cams <- merge(ncl_cams, ca_vegmet, by <- "site", all.x <- TRUE)
ncl_cams$date <- dmy(ncl_cams$date)
ncl_cams$week <- as.numeric(format(ncl_cams$date, "%U")) + 1  # week number
#ncl_cams$week <- as.Date(ncl_cams$date)
ncl_cams$site <- paste(ncl_cams$site, ncl_cams$week, sep <- "_")

ncl_aud <- merge(ncl_aud, ca_vegmet, by <- "site", all.x <- TRUE)
ncl_aud$date <- dmy(ncl_aud$date)
ncl_aud$week <- as.numeric(format(ncl_aud$date, "%U")) + 1  # week number
#ncl_aud$week <- ncl_aud$date
ncl_aud$site <- paste(ncl_aud$site, ncl_aud$week, sep <- "_")

ncl_camaud_week <- rbind(ncl_cams, ncl_aud)

ncl_week_vegmet <- ncl_camaud_week[, c("site", "type", "habitat", "location")]
names(ncl_week_vegmet)[names(ncl_week_vegmet) == "site"] <- "sample"
ncl_week_vegmet <- unique(ncl_week_vegmet)
# add in edna metadata
ncl_week_vegmet <- rbind(ncl_week_vegmet, nclvegmet_fil)

# Remove the scientific_name column (or any other unnecessary columns)
ncl_camaud_week <- ncl_camaud_week[, c("site", "species")]
unique(ncl_camaud_week$species) # should be 64

# method matrix by week
names(ncl_camaud_week)[names(ncl_camaud_week) == "site"] <- "sample"
ncl_camaud_week <- as.data.frame(ncl_camaud_week)

## make abundance matrix for bc nmds plots by habitat
ncl_camaud_week_ab <- as.data.frame(table(ncl_camaud_week$sample, ncl_camaud_week$species))
colnames(ncl_camaud_week_ab) <- c("sample", "species", "abundance")
colnames(ncl_prc_long) <- c("sample", "species", "abundance")
ncl_all_week_ab <- rbind(ncl_camaud_week_ab, ncl_prc_long)
ncl_all_week_ab_wide <- reshape(ncl_all_week_ab, idvar <- "sample", timevar <- "species", direction <- "wide")
colnames(ncl_all_week_ab_wide) <- gsub("^abundance.", "", colnames(ncl_all_week_ab_wide))
ncl_all_week_ab_wide[is.na(ncl_all_week_ab_wide)] <- 0
rownames(ncl_all_week_ab_wide) <- ncl_all_week_ab_wide$sample
ncl_all_week_ab_wide <- ncl_all_week_ab_wide[, -which(names(ncl_all_week_ab_wide) == "sample")]

write.csv(ncl_all_week_ab_wide,'ncl_abundance_matrix_week.csv')

# 

ncl_all_week <- rbind(ncl_camaud_week, ncl_long)

ncl_all_pa_week <- with(ncl_all_week, table(sample, species))

ncl_all_pa_week <- as.data.frame.matrix(ncl_all_pa_week) # should have 83 variables

# presence/absence:
ncl_all_pa_week[ncl_all_pa_week > 0] <- 1
# Remove rows with a sum of zero
empty_samples_ca <- ncl_all_pa_week[rowSums(ncl_all_pa_week) == 0, ]

ncl_all_pa_week <- ncl_all_pa_week[rowSums(ncl_all_pa_week) > 0, ]

ncl_all_pa_week$Columba_oenas <- ifelse(ncl_all_pa_week$Columba_oenas == 1 | ncl_all_pa_week$Columba == 1, 1, 0)
ncl_all_pa_week[["Columba_oenas"]] <- NULL

ncl_all_pa_week$Columba_palumbus <- ifelse(ncl_all_pa_week$Columba_palumbus == 1 | ncl_all_pa_week$Columba == 1, 1, 0)
ncl_all_pa_week[["Columba_palumbus"]] <- NULL

ncl_all_pa_week$Parus_major <- ifelse(ncl_all_pa_week$Parus_major == 1 | ncl_all_pa_week$Paridae == 1, 1, 0)
ncl_all_pa_week[["Paridae"]] <- NULL

ncl_all_pa_week$Anatidae <- ifelse(ncl_all_pa_week$Anatidae == 1 | ncl_all_pa_week$Anas_crecca == 1, 1, 0)
ncl_all_pa_week[["Anas_crecca"]] <- NULL

ncl_all_pa_week[["Passeriformes"]] <- NULL # too vague for current analysis

ncl_all_pa_week <- ncl_all_pa_week[rowSums(ncl_all_pa_week) > 0, ]

ncl_all_pa_week$Phasianidae <- pmax(ncl_all_pa_week$Phasianidae, 
                                    ncl_all_pa_week$Phasianus_colchicus, 
                                    na.rm <- TRUE)
ncl_all_pa_week <- ncl_all_pa_week[, colnames(ncl_all_pa_week) != "Phasianus_colchicus"]

write.csv(ncl_all_pa_week, 'all_week_pa.csv')

# Keep only the rows in the metadata that have corresponding samples in the presence-absence matrix
ncl_week_vegmet <- ncl_week_vegmet[ncl_week_vegmet$sample %in% rownames(ncl_all_pa_week), ]
write.csv(ncl_week_vegmet, "all_week_meta.csv")

# ---- Load cam and acoustic data ----
ncl_long <- fread('ncl_long.csv')
ncl_cams <- fread("ncl_cams.csv")
ncl_aud  <- fread("ncl_acoust.csv")

ncl_cams[, species := gsub(" ", "_", species)]
ncl_aud[, species := gsub(" ", "_", species)]

ncl_cams[, time := ifelse(grepl("^\\d{2}:\\d{2}$", time), paste0(time, ":00"), time)]
ncl_aud[, time := ifelse(grepl("^\\d{2}:\\d{2}$", time),  paste0(time, ":00"), time)]

ncl_cams[, datetime := as.POSIXct(paste(dmy(date), time), format <- "%Y-%m-%d %H:%M:%S")]
ncl_aud[, datetime := as.POSIXct(paste(dmy(date), time), format <- "%Y-%m-%d %H:%M:%S")]

ncl_cams[, type := "Cam"]
ncl_aud[, type := "Aud"]

ncl_combined <- rbindlist(list(ncl_cams, ncl_aud), use.names <- TRUE)

# ± 7 day window
window_days <- 7

# Convert to data.table
setDT(ncl_combined)

# Define rolling intervals
x <- ncl_combined[, .(
  focal_species <- species,
  site,
  datetime_start <- datetime - days(window_days),
  datetime_end   <- datetime + days(window_days),
  focal_time     <- datetime
)]

y <- ncl_combined[, .(
  species,
  site,
  datetime_start <- datetime,
  datetime_end   <- datetime
)]

# Set keys and do foverlaps
setkey(x, site, datetime_start, datetime_end)
setkey(y, site, datetime_start, datetime_end)

co_occurrences <- foverlaps(y, x, type <- "any", nomatch <- 0L)

# Assign sample ID as site + focal_time (date only)
co_occurrences[, sample := paste0(site, "_", format(focal_time, "%Y-%m-%d"))]

# Make long-format presence table
pa_camaud_long <- unique(co_occurrences[, .(sample, species)])
pa_camaud_long[, presence := 1]

# eDNA long presence
ncl_long_pa <- ncl_long[
  , .(presence <- as.integer(sum(reads > 0) > 0)), by <- .(sample, species)
][presence == 1]

# Combine both long-format P/A
combined_long <- rbindlist(list(ncl_long_pa, pa_camaud_long), use.names <- TRUE, fill <- TRUE)
combined_long[is.na(presence), presence := 1]  # Just in case any missing values snuck in

# Remove duplicates (in case species detected by both methods)
combined_long <- unique(combined_long[, .(sample, species, presence)])

# Wide-format PA matrix
combined_pa <- dcast(combined_long, sample ~ species, value.var <- "presence", fill <- 0)

#--- Taxonomic consolidation for combined_pa (same as ncl_all_pa_week) ---
consolidate_taxa_wide <- function(df) {
  # Consolidate pigeon taxa
  if (all(c("Columba", "Columba_oenas") %in% names(df))) {
    df$Columba <- pmax(df$Columba, df$Columba_oenas, na.rm <- TRUE)
    df$Columba_oenas <- NULL
  }
  if (all(c("Columba", "Columba_palumbus") %in% names(df))) {
    df$Columba <- pmax(df$Columba, df$Columba_palumbus, na.rm <- TRUE)
    df$Columba_palumbus <- NULL
  }
  
  # Consolidate Paridae to Parus_major
  if (all(c("Paridae", "Parus_major") %in% names(df))) {
    df$Parus_major <- pmax(df$Parus_major, df$Paridae, na.rm <- TRUE)
    df$Paridae <- NULL
  }
  
  # Consolidate Phasianidae to Phasianus_colchicus
  if (all(c("Phasianidae", "Phasianus_colchicus") %in% names(df))) {
    df$Phasianus_colchicus <- pmax(df$Phasianidae, df$Phasianus_colchicus, na.rm <- TRUE)
    df$Phasianidae <- NULL
  }
  
  # Consolidate Anas_crecca to Anatidae
  if (all(c("Anatidae", "Anas_crecca") %in% names(df))) {
    df$Anatidae <- pmax(df$Anas_crecca, df$Anatidae, na.rm <- TRUE)
    df$Anas_crecca <- NULL
  }
  
  # Drop vague higher-order group
  df$Passeriformes <- NULL
  
  return(df)
}

combined_pa <- consolidate_taxa_wide(combined_pa)


# Output
write.csv(combined_pa, "cooccurrence_pa.csv", row.names <- FALSE)
