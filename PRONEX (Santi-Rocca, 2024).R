########################################################################
#                                PRONEX                                #
########################################################################

##################### ORIGINAL DATASET #################################

### To be applied to an original datasel in .csv format with 4 variables: 
### "FID", "species_name", "taxonomy", "num_reads"
### 'FID' = The series of unique reads come out from sequencing, and for each of them:
### 'species_name' = each species name attributed to the FID with the highest confidence level (superior to 98% in any case is recommended)
###  NOTA: Unique reads will be repeated in the first column as many times as there are species names attributed with the same highest level of confidence)
### 'taxonomy' : Full taxonomic description from phylum to corresponding taxon level
### 'num_reads' : number of ocurrences of each unique FID

#### User Input: Thresholds for Cluster Importance and Taxon Abundance ####

## Define here the thresholds:
# Fraction of total reads for cluster fusion
FID_importance_fraction <- 1 / 10000  # Default value
# Fraction of total reads for permanent exclusion of FIDs
exclusion_threshold <- 1 / 200000  # Default value

########################################################################
#          Cite the initial script (Santi-Rocca et al. 2024)          #
########################################################################

# Load necessary library
library(dplyr)

#### CUSTOMIZE THIS SECTION ####
# Set working directory (replace with your own directory)
# Example: setwd("C:/path/to/your/working/directory")
setwd("YOUR_WORKING_DIRECTORY_HERE")

# Specify the input CSV file name (replace with your file name)
input_file <- "YOUR_INPUT_FILE_HERE.csv"

# Specify the output CSV file name
output_file <- "Results.csv"

#################################

# Read the data file
data <- read.csv(input_file, header = TRUE)

# Expected column names (update as necessary to match your dataset)
expected_columns <- c("FID", "species_name", "taxonomy", "num_reads")

# Check if the column names are as expected
if (!all(expected_columns %in% names(data))) {
  stop(
    "Column names do not match the expected format. Expected columns are: ",
    paste(expected_columns, collapse = ", ")
  )
}

# Calculate total number of reads (without duplicates for FIDs)
total_reads <- data %>%
  distinct(FID, .keep_all = TRUE) %>%
  summarise(total_reads = sum(num_reads, na.rm = TRUE)) %>%
  pull(total_reads)

# Calculate thresholds based on total reads
FID_importance_threshold <- total_reads * FID_importance_fraction
exclusion_threshold_absolute <- total_reads * exclusion_threshold

########################################################################
# Step 1: Exclude FIDs permanently based on exclusion threshold       #
########################################################################

data <- data %>%
  filter(num_reads >= exclusion_threshold_absolute)

########################################################################
# Initialize columns for clustering and thresholds                   #
########################################################################

# Add columns for cluster importance and taxon abundance
data <- data %>%
  mutate(
    initial_cluster = NA,
    cluster_importance = ifelse(num_reads >= FID_importance_threshold, 1, 0),
    taxon_abundance = ifelse(num_reads >= FID_importance_threshold, 1, 0)
  )

########################################################################
# PHASE 1: Initialize clustering and fuse clusters                   #
########################################################################

# Initialize clustering based on FID occurrences
data$initial_cluster[1] <- 1

for (i in 2:nrow(data)) {
  if (sum(data$FID[1:(i - 1)] == data$FID[i]) == 0) {
    data$initial_cluster[i] <- max(data$initial_cluster[1:(i - 1)]) + 1
  } else {
    data$initial_cluster[i] <- min(data$initial_cluster[1:(i - 1)][data$FID[1:(i - 1)] == data$FID[i]])
  }
}

# Define a function to fuse clusters by species overlap
fuse_clusters_by_species <- function(data) {
  data <- data %>%
    group_by(species_name) %>%
    mutate(
      lowest_cluster = ifelse(cluster_importance == 1, min(initial_cluster, na.rm = TRUE), initial_cluster)
    ) %>%
    ungroup() %>%
    group_by(initial_cluster) %>%
    mutate(initial_cluster = min(lowest_cluster, na.rm = TRUE)) %>%
    ungroup() %>%
    select(-lowest_cluster)
  
  return(data)
}

# Iteratively fuse clusters based on species overlap and cluster importance
iteration <- 1
clusters_changed <- TRUE

while (clusters_changed) {
  new_data <- fuse_clusters_by_species(data)
  clusters_changed <- !all(new_data$initial_cluster == data$initial_cluster)
  if (clusters_changed) {
    iteration <- iteration + 1
    data <- new_data
    colname <- paste0("iteration_", iteration)
    data[[colname]] <- data$initial_cluster
  }
}

########################################################################
# Phase 2: Subdivide and refine clusters                             #
########################################################################

# Add columns for further refinement
data <- data %>%
  group_by(species_name) %>%
  mutate(Important_species = ifelse(any(cluster_importance == 1), 1, 0)) %>%
  ungroup() %>%
  group_by(FID) %>%
  mutate(Only_not_important = ifelse(all(Important_species == 0), 1, 0)) %>%
  ungroup() %>%
  mutate(To_subdivide = Important_species + Only_not_important) %>%
  group_by(species_name) %>%
  mutate(Phase_2 = ifelse(To_subdivide != 0, min(initial_cluster, na.rm = TRUE), NA)) %>%
  ungroup() %>%
  group_by(FID) %>%
  mutate(Division_factor = sum(To_subdivide > 0, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Division_factor = ifelse(is.na(Phase_2), 0, Division_factor))

########################################################################
# Save the final dataset                                              #
########################################################################

# Save the resulting dataset to a CSV file
write.csv(data, output_file, row.names = FALSE)

# Inform the user
cat("The script has completed successfully. Results saved to", output_file, "\n")
