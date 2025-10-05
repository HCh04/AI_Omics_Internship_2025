#####################################################################
# HasmikChilingaryan_ClassIII_Preprocess.R
# Module III — Class 3B: Preprocessing & Normalisation (E-MEXP-1277)
# Workflow: import raw text files (Agilent), QC, normalize, filter,
# produce boxplot & PCA images (for submission), save results.
#####################################################################

# Check working directory

getwd()
list.files()

# Library Installation and Load
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c("limma", "arrayQualityMetrics", "Biobase", "genefilter"), ask = FALSE, update = TRUE)

# CRAN packages
install.packages(c("matrixStats", "ggplot2", "dplyr"), dependencies = TRUE, quiet = TRUE)

# Load packages
library(limma)
library(arrayQualityMetrics)
library(Biobase)
library(genefilter)
library(matrixStats)
library(ggplot2)
library(dplyr)

# Paths
module_root <- "Module_III_E-MEXP-1277"
raw_dir <- file.path(module_root, "raw_data")
results_dir <- file.path(module_root, "results")
figures_dir <- file.path(module_root, "figures")
scripts_dir <- file.path(module_root, "scripts")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

cat("raw_dir:", raw_dir, "\n")
cat("Files in raw_dir:\n") ; print(list.files(raw_dir))

# Read SDRF (metadata) and inspect group column

sdrf_files <- list.files(raw_dir, pattern = "sdrf", ignore.case = TRUE, full.names = TRUE)
if(length(sdrf_files) == 0) stop("No SDRF file found in raw_data. Upload E-MEXP-1277.sdrf.txt")

sdrf <- read.delim(sdrf_files[1], stringsAsFactors = FALSE, check.names = FALSE)
cat("SDRF loaded/ Rows:", nrow(sdrf), "Columns:", ncol(sdrf), "\n")
cat("SDRF column names:\n"); print(names(sdrf))
cat("Preview SDRF:\n"); print(head(sdrf))

# Find a likely group column automatically
candidate_cols <- names(sdrf)[grep("expos|dose|exposed|control|status|case|group|condition|characteristics", names(sdrf), ignore.case = TRUE)]
cat("Candidate SDRF columns for group info:\n"); print(candidate_cols)

# Show unique values of all candidates to help choose correct column
if(length(candidate_cols) > 0){
  cat("Unique values in candidate SDRF columns:\n")
  for (col in candidate_cols) {
    cat("\nColumn:", col, "\n")
    print(unique(sdrf[[col]]))
  }
  # Pick the first candidate by default; user can change
  chosen_col <- candidate_cols[1]
  cat("\nUsing SDRF column for group labeling:", chosen_col, "\n")
} else {
  stop("No obvious group column found automatically. Inspect SDRF column names and choose one manually.")
}

# Read raw microarray files (Agilent text files) (robust ver.)
raw_txt_files <- list.files(raw_dir, pattern = "\\.txt$", full.names = TRUE)
# Exclude SDRF file
raw_txt_files <- raw_txt_files[!grepl("sdrf", raw_txt_files, ignore.case = TRUE)]
if(length(raw_txt_files) == 0) stop("No Agilent raw text files found in raw_data.")

raw_basenames <- basename(raw_txt_files)

# read.maimages for Agilent single-channel
ma <- read.maimages(
  files = raw_basenames,
  path = raw_dir,
  source = "agilent",
  green.only = TRUE,
  columns = list(
    G  = "Feature Extraction Software:gProcessedSignal",
    Gb = "Feature Extraction Software:gBGMedianSignal"
  )
)

cat("Expression matrix (raw) dimensions (probes x arrays):", dim(ma$E), "\n")
exprs_raw <- ma$E

# QC before Normalisation
# Median z-score
sample_medians_pre <- apply(exprs_raw, 2, median, na.rm = TRUE)
med_z_pre <- scale(sample_medians_pre)
outliers_pre <- names(sample_medians_pre)[which(abs(med_z_pre)>2)]
cat("Median z-score outliers (pre-normalization):\n"); print(outliers_pre)

# Distance-based z-score
dist_mat_pre <- as.matrix(dist(t(exprs_raw)))
avg_dist_pre <- apply(dist_mat_pre,1,mean)
dist_z_pre <- scale(avg_dist_pre)
outliers_dist_pre <- names(avg_dist_pre)[which(abs(dist_z_pre)>2)]
cat("Distance-based outliers (pre-normalization):\n"); print(outliers_dist_pre)

# Boxplot
pdf(file.path(figures_dir,"boxplot_raw_data.pdf"), width=10, height=6)
boxplot(exprs_raw, main="Raw intensity distribution (per array)", las=2, cex.axis=0.6)
dev.off()

# arrayQualityMetrics report
raw_eset <- ExpressionSet(assayData=exprs_raw)
arrayQualityMetrics(expressionset=raw_eset,
                    outdir=file.path(results_dir,"QC_Raw"),
                    force=TRUE, do.logtransform=TRUE)

# Normalisation (single-channel Agilent)
log_expr <- log2(exprs_raw + 1)
norm_mat <- normalizeBetweenArrays(log_expr, method = "quantile")

cat("Normalization complete. Dimensions (probes x arrays):", dim(norm_mat), "\n")

# QC after Normalisation
sample_medians_post <- apply(norm_mat,2,median,na.rm=TRUE)
med_z_post <- scale(sample_medians_post)
outliers_post <- names(sample_medians_post)[which(abs(med_z_post)>2)]

dist_mat_post <- as.matrix(dist(t(norm_mat)))
avg_dist_post <- apply(dist_mat_post,1,mean)
dist_z_post <- scale(avg_dist_post)
outliers_dist_post <- names(avg_dist_post)[which(abs(dist_z_post)>2)]

pdf(file.path(figures_dir,"boxplot_normalized_data.pdf"), width=10, height=6)
boxplot(norm_mat, main="Normalized distribution", las=2, cex.axis=0.6)
dev.off()

norm_eset <- ExpressionSet(assayData=norm_mat)
arrayQualityMetrics(expressionset=norm_eset,
                    outdir=file.path(results_dir,"QC_Normalized"),
                    force=TRUE, do.logtransform=FALSE)

# Filtering Low-Intensity probes
row_medians <- rowMedians(as.matrix(norm_mat))
png(file.path(figures_dir, "median_intensity_distribution.png"), width = 800, height = 600)
hist(row_medians, breaks = 100, main = "Median intensity distribution (norm)", xlab = "Median intensity")
abline(v = 3.5, col = "red", lwd = 2)
dev.off()

threshold <- 3.5
keep_probe_idx <- which(row_medians > threshold)
filtered_mat <- norm_mat[keep_probe_idx, , drop = FALSE]
cat("Probes before filtering:", nrow(norm_mat), "\n")
cat("Probes after filtering :", nrow(filtered_mat), "\n")
write.csv(filtered_mat, file = file.path(results_dir, "filtered_expression.csv"))

# PCA plot (post-filtering)
pca <- prcomp(t(filtered_mat), scale. = TRUE)
pca_df <- as.data.frame(pca$x[,1:2])
pca_df$Sample <- colnames(filtered_mat)

pdf(file.path(figures_dir, "PCA_after_normalization.pdf"), width = 7, height = 6)
plot(pca_df$PC1, pca_df$PC2, xlab = "PC1", ylab = "PC2", main = "PCA after normalization and filtering", pch = 19)
text(pca_df$PC1, pca_df$PC2, labels = pca_df$Sample, cex = 0.6, pos = 3)
dev.off()

# Phenotype Groups
if(exists("chosen_col") && chosen_col %in% names(sdrf)){
  raw_group_values <- sdrf[[chosen_col]]
  group_labels <- ifelse(grepl("expos|radi|xray|x-ray|worker|exposed", raw_group_values, ignore.case = TRUE),
                         "Exposed", "Control")
} else {
  nS <- ncol(norm_mat)
  group_labels <- c(rep("Exposed", floor(nS/2)), rep("Control", nS - floor(nS/2)))
}

groups <- factor(group_labels, levels = c("Control","Exposed"))
metadata <- data.frame(Sample = colnames(norm_mat), Group = groups, stringsAsFactors = FALSE)
write.csv(metadata, file = file.path(results_dir, "phenotype_metadata.csv"), row.names = FALSE)

# Summary for form submission
total_samples <- ncol(norm_mat)
disease_count <- sum(groups=="Exposed")
control_count <- sum(groups=="Control")
probes_before <- nrow(norm_mat)
probes_after <- nrow(filtered_mat)
outliers_before_count <- length(unique(c(outliers_pre,outliers_dist_pre)))
outliers_after_count <- length(unique(c(outliers_post,outliers_dist_post)))

cat("\n--- FORM VALUES ---\n")
cat("Accession: E-MEXP-1277\nTotal samples:", total_samples,
    "\nDisease (Exposed):", disease_count,
    "\nControl:", control_count,
    "\nProbes before filtering:", probes_before,
    "\nProbes after filtering:", probes_after,
    "\nOutliers before normalization:", outliers_before_count,
    "\nOutliers after normalization:", outliers_after_count,"\n")

# Save workspace & summary CSV
save.image(file=file.path(results_dir,"ClassIII_Assignment.RData"))
summary_df <- data.frame(
  accession="E-MEXP-1277",
  total_samples=total_samples,
  disease_samples=disease_count,
  control_samples=control_count,
  probes_before=probes_before,
  probes_after=probes_after,
  outliers_before=outliers_before_count,
  outliers_after=outliers_after_count,
  stringsAsFactors=FALSE
)
write.csv(summary_df, file=file.path(results_dir,"E-MEXP-1277_summary_for_submission.csv"), row.names=FALSE)

cat("FINISHED: check 'results/' and 'figures/' for outputs and QC folders.\n")