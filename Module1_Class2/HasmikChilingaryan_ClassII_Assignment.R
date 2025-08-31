# ===================================================================
#             AI and Omics Research Internship (2025)
#             Module I: Basic R Functions (Class 2)
#            #### notes made for me to follow easier ####
# ===================================================================

# Assignment 2
# --------------------------
# In this assignment you will work with the results of differential gene expression (DGE) analysis. 
#The analysis produces two key measures for each gene:

# log2FoldChange (log2FC): 
# Indicates the magnitude and direction of change in gene expression. 
# Positive values suggest higher expression(upregulated gene) in the experimental condition compared to control. 
# Negative values suggest lower expression (downregulated gene). 
# The absolute value reflects the strength of the change.

# Adjusted p-value (padj): 
# Represents the statistical significance of the observed difference, corrected for multiple testing. 
# A smaller value indicates stronger evidence that the observed difference is not due to chance.

# Write a function classify_gene() 

# that takes:
#   - logFC (log2FoldChange)
#   - padj  (adjusted p-value)

# and returns:
#   - "Upregulated" if log2FC > 1 and padj < 0.05
#   - "Downregulated" if log2FC < -1 and padj < 0.05
#   - "Not_Significant" otherwise


# Then:
#   - Apply it in a for-loop to process both datasets (DEGs_data_1.csv, DEGs_data_2.csv)
#   - Replace missing padj values with 1
#   - Add a new column 'status'
#   - Save processed files into Results folder
#   - Print summary counts of significant, upregulated, and downregulated genes
#   - Use table() for summaries

# Data Availability
# The input files are available in the GitHub repository:
#      DEGs_Data_1.csv
#      DEGs_Data_2.csv

# Each file contains three columns: 
# Gene_Id	
# padj	
# logFC


# Setup input/output

input_dir <- "raw_data" #folder where CSVs are uploaded
output_dir <- "results" #folder to save outputs

# Make sure the results folder exists
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# List files to process
files_to_process <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

# Function to classify_gene() takes logFC and padj 
# and returns gene status

# Define classify_gene() function
classify_gene <- function(logFC, padj) {
  if (is.na(padj) & padj < 0.05 & logFC > 1) {
    return("Upregulated")
  } else if(is.na(padj) & padj < 0.05 & logFC < -1) {
    return("Downregulated")
  } else {
    return("Not_Significant")
  }
}

# Empty list to store processed results
result_list <- list()

# Process each file in a for-loop
for(file_name in files_to_process) {
  cat("\nProcessing:", file_name, "\n") #print progress
  
  # Build the full file path for input
  input_file_path <- file.path(input_dir, file_name)
  
  # Import DEGs dataset
  degs_data <- read.csv(input_file_path, header = TRUE)
  cat("File imported. Columns detected:\n")
  print(names(degs_data))
  
  # Replace missing padj values with 1
  degs_data$padj[is.na(degs_data$padj)] <- 1
  
  # Apply classification fuction across all rows
  degs_data$status <- mapply(classify_gene,
                             logFC = degs_data$logFC,
                             padj = degs_data$padj)
  
  # Save results in R list
  result_list[[file_name]] <- degs_data
  
  # Save processed file into results folder
  output_file_path <- file.path(output_dir, paste0("Processed_", file_name))
  write.csv(degs_data, output_file_path, row.names = FALSE)
  
  # Print summary counts for this dataset
  cat("Summary of gene status:\n")
  print(table(degs_data$status))
}

result_list[["DEGs_data_1.csv"]]
result_list[["DEGs_data_2.csv"]]

save.image("HasmikChilingaryan_ClassIIb_Assignment.RData")

