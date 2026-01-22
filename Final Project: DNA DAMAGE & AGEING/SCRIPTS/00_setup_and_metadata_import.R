############################################################
# 00_setup_and_metadata_import.R
# Project: DNA Damage & Ageing Final Project
# Step: Set up directories & import metadata for GSE206778
############################################################

# Load packages
library(tidyverse)

# Define project directory
proj_dir <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA"

# Define data directories
raw_dir  <- file.path(proj_dir, "data_raw")
qc_dir   <- file.path(proj_dir, "data_qc")
res_dir  <- file.path(proj_dir, "results")
script_dir <- file.path(proj_dir, "scripts")

# Create folders if not existing
dir.create(raw_dir, showWarnings = FALSE)
dir.create(qc_dir, showWarnings = FALSE)
dir.create(res_dir, showWarnings = FALSE)
dir.create(script_dir, showWarnings = FALSE)

# Path to SRA metadata for GSE206778
meta_file <- file.path(raw_dir, "GSE206778_SRA_metadata.csv")

# Load metadata
meta <- read_csv(meta_file)

# Print summary
cat("Metadata imported successfully.\n")
cat("Number of samples: ", nrow(meta), "\n")
cat("Columns available:\n")
print(colnames(meta))

# Save a cleaned version (optional)
meta_clean <- meta %>%
  rename(sample = Run) %>%  # Ensures a standard name for sample IDs
  mutate(sample = as.character(sample))

write_csv(meta_clean, file.path(raw_dir, "GSE206778_metadata_clean.csv"))

cat("Cleaned metadata saved to: data_raw/GSE206778_metadata_clean.csv\n")
