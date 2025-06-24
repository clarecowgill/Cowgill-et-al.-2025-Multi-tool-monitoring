README
================
2025-06-23

# Multi-tool monitoring: integrating eDNA, acoustics, and cameras to monitor taxonomic, functional and trophic diversity

## Table of Contents

- [Overview](#overview)
- [Dataset Contents](#dataset-contents)
  - [CSV Data Files](#csv-data-files)
  - [Media Files](#media-files)
  - [DNA Sequences](#dna-sequences)
  - [Data Analysis Scripts](#data-analysis-scripts)
- [Data Collection and Methods](#data-collection-and-methods)
- [Usage Notes](#usage-notes)
- [Citation](#citation)
- [Contact](#contact)

## Overview

This dataset supports a multi-method biodiversity assessment combining
environmental DNA (eDNA) metabarcoding with passive acoustic and camera
trapping methods. Data was collected at a rewilding site in Scotland to
compare taxonomic richness, functional diversity, and trophic
co-occurrence networks detected by these complementary methods.

## Dataset Contents

### CSV Data Files

- `ncl_matrix_raw.csv`  
  Raw eDNA sequence-based detections from water, soil, tree rolling, and
  faecal samples.
  
  _Note on sample naming: the first letter of each sample corresponds to substrate type (W=Water, T=Tree roller, S=Soil, P=Scat), the second correpsonds to habitat (C=Confier, B=Broadleaf, F=Felled, G=Bog) with sample number following. For soil and scat samples there is a suffix of either 1 or 2 denoting that these are lab replicates of the same field sample. Blanks are labelled as: NEG (PCR negative), POS (PCR positive), 'EB' (e.g. TEB1) are extraction blanks, 'BK' is a field blank (e.g. TBK1), FB (e.g. TFB1) denotes a filter blank when different distilled water was used._

- `blank_associations.csv`  
  Field, extraction, and PCR blanks and the samples they are associated
  with- used to create read thresholds to prevent fale positives.

- `all_week_pa.csv`  
  Presence-absence matrix summarising detections across all methods.

- `all_week_meta.csv`  
  Metadata for all eDNA samples and each week of acoustic and camera
  monitoring.

- `ncl_acoust.csv`  
  Acoustic data from AudioMoth devices, processed with BirdNet and BTO
  Acoustic Pipeline.

- `ncl_cams.csv`  
  Camera trap detection data.

- `nodes.csv`  
  Functional trait data for detected species (body mass, trophic guild,
  lifestyle). Trait data collected from multiple open source databases:
  Pantheria, EUBirds, Amphibio, Avonet.

- `cooccurrence_pa.csv`  
  Presence-absence matrix across all methods with dynamic week-long
  co-occurrence data. This dataframe is used for the network analysis to
  enable dynamic week-long co-occurences (i.e. co-occurrence of ± 7 days
  for each identification as opposed to strict week-long units).

- `cooccurrence_meta.csv`  
  Metadata corresponding to the co-occurrence presence-absence matrix.

- `trophic_links.csv`  
  Database of trophic links used to construct trophic metawebs. Trophic
  link information primarily gathered from GLOBI database, with
  additional information collected to supplement gaps. References given
  for these additions.

### Media Files

- Zipped acoustic recordings from 7 AudioMoth devices and  camera trap photos from 15 cameras at 8 locations are available on the Zenodo repository for this publication

### DNA Sequences

The Zenodo repository for this paper contains the raw DNA sequence data generated from 10 sequencing libraries using the Illumina MiSeq platform. 

- The full dataset is provided as a single compressed file: `DNA Sequences.zip`.

- When unzipped, the archive reveals 10 folders named `Lib1` through `Lib10`, each containing raw paired-end sequencing files in compressed FASTQ format (`.fastq.gz`).

- Each FASTQ file includes sequence and quality score information and follows Illumina’s standard naming conventions (e.g., `Lib1_S1_L001_R1_001.fastq.gz` and      `Lib1_S1_L001_R2_001.fastq.gz`).

- The FASTQ files are compressed using BGZF (Blocked GNU Zip Format), compatible with standard gzip tools and bioinformatics software.

### Reference Database

Database used for taxonomic identification of vertebrates in fasta format (`12S_verts.fasta`) and txt file mapping sequences to their taxonomic IDs from Genbank (`12S_verts_tax_map.txt`).

### Data Analysis Scripts

  ------------------------------------------------------------------------

`initial_cleanup.R`  

Data cleaning and preparation pipeline: filters contamination, standardises taxonomy,
creates presence-absence matrices and metadata.

##### Input Files
  
  - `ncl_raw.tsv` — Raw OTU table for eDNA samples  
  - `blank_association_all.csv` — Blank/sample associations  
  - `edna_meta.csv` — Metadata for eDNA samples  
  - `ncl_cams.csv` — Camera trap detection data  
  - `ncl_acoust.csv` — Acoustic detection data  
  - `all_meta.csv` — Combined metadata across methods
  
  
#### Workflow Summary
  
  **1. Clean and Format Raw eDNA Data**
  
  - Removes unassigned taxa and non-NCL columns
  - Drops OTUs with zero reads  
  - **Output:** `ncl_matrix_raw.csv`
  
  **2. Filter Contamination (LOD/LOQ)**
  
  - Merges blank associations and metadata
  - Calculates limit of detection (LOD) and limit of quantification (LOQ)
  - Filters out reads below LOD
  - **Output:** filtered species-by-sample matrix with taxonomy
    (`ncl_filtered.csv`)
  
  **3. Remove Non-target Taxa and Controls**
  
  - Excludes fish, domestic animals, and control samples (e.g. NEG, POS)
  - Applies 0.1% abundance threshold to remove trace reads
  - **Output:** Updated matrix (`ncl_filtered.csv`)
  
  **4. Taxonomic Consolidation**
  
  - Combines overlapping taxa (e.g. genus-level and species-level rows)
  - E.g. merges `"Columba_oenas"` and `"Columba_palumbus"` into
    `"Columba"`
  - **Output:** `ncl_f.csv`
  
  **5. Prepare Matrices for Diversity Analyses**
  
  - Converts read data to:
    - Presence/absence matrix (`ncl_pa.csv`)
    - Proportional abundance matrix (`ncl_prc.csv`)
  - Filters out empty samples
  - Merges with metadata (`nclvegmet_fil.csv`)
  - Converts to long format (`ncl_long.csv`)
  
  **6. Combine eDNA, Camera, and Acoustic Data**
  
  - Standardises species names
  - Merges acoustic and camera trap records
  - Creates multi-method presence-absence matrix (`ncl_all_pa.csv`)
  - Removes redundant family/genus-level entries
  
  **7. Weekly Resolution Data for Temporal Analysis**
  
  - Converts camera and acoustic dates to week-based site IDs
    (e.g. `NCL01_18`)
  - Aggregates abundance and presence data by week
  - Merges with eDNA and metadata
  - **Outputs:**
    - Abundance matrix: `ncl_abundance_matrix_week.csv`
    - Presence matrix: `all_week_pa.csv`
    - Metadata: `all_week_meta.csv`
  
  **8. Co-occurrence Analysis**
  
  - Combines eDNA, camera, and acoustic data into single matrix
  - Includes ±7 day temporal matching for camera/acoustic records
  - Outputs long- and wide-format presence data for co-occurrence
    inference
  - Final matrix: `cooccurrence_pa.csv`

------------------------------------------------------------------------

`ncl_model.R`  
  Generalized linear models assessing drivers of species richness. Comparing different variables and their influence on species richness within each sample. Looks at model AIC scores to determine the best fit and creates diagnostic plots.

------------------------------------------------------------------------

`ncl_vegan.R`  
  
  Alpha and beta diversity analyses using the vegan R package. This script calculates diversity metrics (species richness, jaccard dissimilarity), generates plots (boxplots for alpha diversity, NMDS for community diversity), and runs statistical tests to compare community composition across sample types and detection methods.
    
  **Ouputs**
  
  - Plots: Boxplots and NMDS ordination plots

  - Statistical summaries: results from Dunn’s test, ANOVA, PERMANOVA.
  
------------------------------------------------------------------------

`specaccums.R`  
  Species accumulation and extrapolation analysis with iNEXT package.

  **Ouputs**
  
  - Plots: Species accumulation curves, extrapolated to faciliate comparisons:
    a) DNA samples, along with an 'all' sample types (pooling all methods) line
    b) Passive methods, with taxonomic richness against the number of days (pooled across multiple sites)

------------------------------------------------------------------------

`habitat_discrimination.R`  
  Understanding community differences between habitats, looking individually at each method. NMDS plots looking at Jaccard dissimilarity and testing differences via pairwise PERMANOVA.
    
  **Ouputs**
  
  - Plots: NMDS plots of data split by method, looking at the differences in community compared to the different habitats samples were collected from

  - Statistical analyses: Pairwise PERMANOVAs for each habitat pairwise comparison
  
------------------------------------------------------------------------

`trait_analysis.R`  
  Functional diversity analyses (taxonomic distinctiveness, Rao’s
  quadratic entropy, functional redundancy). Also visualises traits in space relative to other samples with PCoA.
    
  **Ouputs**
  
  - Plots:
    a) Boxplots with points for each eDNA sample / week of passive sampling. Comparing taxonomic distinctiveness, Rao's quadratic entropy, Functional redundancy
    b) PCoA plot of traits and points for sampling units

  - Statistical analyses: Dunn's tests to compare functional diversity metrics
    
------------------------------------------------------------------------

`networks.R`  
  Construction of trophic metawebs by combining co-occurrence and
  trophic link data. Calculates a number of commonly used network indices for all possible combinations of methods (linkage density, connectance, number of links, modularity, robustness and mean trophic level)
      
  **Ouputs**
  
- Plots:
    - Visualise networks for all the detected trophic links, with nodes set in specific positions
    - Bar plots of network metrics

 <small>_Note: to create different network visualisation, it is required that the colour palette highlighting the focal methods is specified, and the file name to save the png is changed for each network._</small>

------------------------------------------------------------------------

## Data Collection and Methods

Sampling was conducted in June 2023 across four microhabitats (broadleaf woodland, conifer woodland, felled and bog) at the
Natural Cqpital Lab rewilding site in Scotland. eDNA samples were collected from water, soil, tree
rolling, and faecal sources and analyzed by metabarcoding. Passive
acoustic monitoring employed AudioMoth recorders, and Browning BTC-8E-H4
camera traps were used for vertebrate detection. Data were integrated to
assess taxonomic richness, functional diversity, and trophic network
completeness.

## Usage Notes

- The files are structured to facilitate joint analysis of eDNA, acoustic, and camera trap data.
- For detailed file descriptions and workflows, please refer to the analysis scripts and inline comments.
- The typical workflow begins with running `initial_cleanup.R`. Following this, any of the scripts can be run with the subsequent outputs.
- Users should ensure required R packages are installed.
- Zenodo repository DOI: 10.5281/zenodo.15723455

