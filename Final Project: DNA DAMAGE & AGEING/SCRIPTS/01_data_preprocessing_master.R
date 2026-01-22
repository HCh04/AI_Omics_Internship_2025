###############################################################
# MASTER SCRIPT — DATA RETRIEVAL, QC, NORMALISATION
# Final Project: DNA Damage & Ageing
###############################################################

# -----------------------------
# Step 0: Load required packages
# -----------------------------
library(tidyverse)
library(DESeq2)
library(pheatmap)

# -----------------------------
# Step 0a: Set up local library
# -----------------------------
# Choose a folder where you have write access
local_lib <- "Rlibs"
dir.create(local_lib, showWarnings = FALSE)

# -----------------------------
# Step 0b: Install Bioconductor manager if missing
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = local_lib)
}

# -----------------------------
# Step 0c: Install Bioconductor packages locally
# -----------------------------
BiocManager::install("zellkonverter", lib = local_lib, ask = FALSE, update = FALSE)
BiocManager::install("SummarizedExperiment", lib = local_lib, ask = FALSE, update = FALSE)

# -----------------------------
# Step 0d: Load Bioconductor packages from local library
# -----------------------------
library(zellkonverter, lib.loc = local_lib)
library(SummarizedExperiment, lib.loc = local_lib)

# -----------------------------
# Step 1: Define file paths
# -----------------------------
data_dir <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/data_raw"

# Files
GSE206778_metadata <- file.path(data_dir, "GSE206778_metadata_clean.csv")
GSE206778_diff     <- file.path(data_dir, "GSE206778_ross_dna_damage_gene_exp.diff")
GSE288213_counts   <- file.path(data_dir, "GSE288213_raw_counts_Csa_WT.csv")
GSE209742_h5ad     <- file.path(data_dir, "GSE209742_Patel_LSK_raw_Submission.h5ad")

# -----------------------------
# Step 2: Load datasets
# -----------------------------
# GSE206778 - differential expression (already processed)
metadata_206778 <- read.csv(GSE206778_metadata, header = TRUE)

# Load data.table
library(data.table)

# Path to your diff file
diff_file <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/data_raw/GSE206778_ross_dna_damage_gene_exp.diff/GSE206778_ross_dna_damage_gene_exp.diff"

# Read only necessary columns
diff_206778 <- fread(diff_file, select = c("gene_id", "`log2(fold_change)`", "p_value", "q_value"))

# Inspect first few rows
head(diff_206778)
dim(diff_206778)

# Read just the first few lines to see headers
diff_header <- fread(diff_file, nrows = 0)
colnames(diff_header)

# GSE288213 - raw counts RNA-seq
GSE288213_counts <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/data_raw/GSE288213_raw_counts_Csa_WT.csv/GSE288213_raw_counts_Csa_WT.csv"

# Load using fread (efficient for large CSV)
counts_288213 <- fread(GSE288213_counts, header = TRUE)

# Convert first column to rownames if needed
counts_288213 <- as.data.frame(counts_288213)
rownames(counts_288213) <- counts_288213[,1]
counts_288213 <- counts_288213[,-1]

# Inspect first few rows
head(counts_288213)
dim(counts_288213)

# GSE209742 - h5ad single-cell or bulk RNA-seq
expr_209742_file <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/data_raw/GSE209742_Patel_LSK_raw_Submission.h5ad/GSE209742_expr_HVG.csv"
expr_209742 <- fread(expr_209742_file, header = TRUE)

# Convert first column to rownames
rownames(expr_209742) <- expr_209742[[1]]
expr_209742 <- expr_209742[,-1]
expr_209742 <- as.data.frame(expr_209742)
head(expr_209742)
dim(expr_209742)

###############################################################
# Step 3: Quality Control & Normalization
###############################################################

# -----------------------------
# 3a. GSE206778 diff data
# -----------------------------
# Already processed DE results; just ensure numeric
diff_206778$`log2(fold_change)` <- as.numeric(diff_206778$`log2(fold_change)`)
diff_206778$p_value <- as.numeric(diff_206778$p_value)
diff_206778$q_value <- as.numeric(diff_206778$q_value)

# Optional: filter significant genes
diff_sig_206778 <- diff_206778 %>%
  filter(!is.na(q_value) & q_value < 0.05)

# -----------------------------
# 3b. GSE288213 counts
# -----------------------------
# Install DESeq2
BiocManager::install("DESeq2")

# Install pheatmap and other CRAN packages
install.packages(c("pheatmap", "tidyverse", "caret", "randomForest", "ggplot2"))
library(DESeq2)

# Convert to DESeq2 object for normalisation
dds_288213 <- DESeqDataSetFromMatrix(
  countData = counts_288213,
  colData = data.frame(row.names = colnames(counts_288213)),
  design = ~1
)

# Filter low counts
keep <- rowSums(counts(dds_288213)) >= 10
dds_288213 <- dds_288213[keep, ]

# Normalise counts
dds_288213 <- estimateSizeFactors(dds_288213)
norm_counts_288213 <- counts(dds_288213, normalized = TRUE)

# -----------------------------
# 3b. CLEAN METADATA — GSE288213
# -----------------------------

library(tidyverse)

# -------------------------------------------
# 1. Define file paths
# -------------------------------------------
proj_dir <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA"
raw_dir  <- file.path(proj_dir, "data_raw")

# This is your actual metadata file path
meta_file <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/data_raw/GSE288213_raw_counts_Csa_WT.csv/SraRunTable.csv"

# -------------------------------------------
# 2. Load metadata
# -------------------------------------------
meta <- read.csv(meta_file, header = TRUE, stringsAsFactors = FALSE)

cat("Metadata loaded.\n")
cat("Dimensions: ", paste(dim(meta), collapse = " x "), "\n")

# -------------------------------------------
# 3. Inspect columns
# -------------------------------------------
print(colnames(meta))

# -------------------------------------------
# 4. Select useful biological columns
# -------------------------------------------
meta_clean <- meta %>%
  select(
    sample  = Run,      
    genotype,
    age = time,           # renamed 'time' to 'age'
    tissue,               # replace with actual tissue column name if different
    BioSample,
    Experiment,
    LibraryLayout,
    Organism
  )

# -------------------------------------------
# 5. Clean fields and convert age to numeric
# -------------------------------------------
meta_clean <- meta_clean %>%
  mutate(
    sample = as.character(sample),
    genotype = trimws(genotype),
    age = as.character(age),
    # Convert age to numeric days
    age_days = case_when(
      grepl("month", age) ~ as.numeric(gsub(" months", "", age)) * 30,
      grepl("day", age) ~ as.numeric(gsub(" days", "", age)),
      TRUE ~ NA_real_
    )
  )

# -------------------------------------------
# 6. Save cleaned metadata
# -------------------------------------------
output_path <- file.path(raw_dir, "GSE288213_metadata_clean.csv")
write.csv(meta_clean, output_path, row.names = FALSE)

cat("Cleaned metadata saved at:\n", output_path, "\n")

library(DESeq2)

# counts_288213 should be raw counts (not normalized)
dds <- DESeqDataSetFromMatrix(
  countData = counts_288213,
  colData = meta,          # make sure meta has rows matching columns of counts
  design = ~ genotype      # or ~ genotype + age_days if you want to include age
)

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res[order(res$padj), ])
head(res)

# Export DEGs
deg_288213_file <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/preprocessed/DEGs_288213.csv"
write.csv(res, deg_288213_file, row.names = TRUE)


# -----------------------------
# 3c. GSE209742 HVG expression
# -----------------------------
# Already selected HVGs, just make sure numeric
expr_209742 <- as.data.frame(lapply(expr_209742, as.numeric))
rownames(expr_209742) <- rownames(expr_209742) # preserve rownames

# Optional: log-transform for comparability
expr_209742_log <- log1p(expr_209742)

# Create preprocessed folder if it doesn't exist
dir.create("U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/preprocessed", showWarnings = FALSE)

# GSE206778 — DE genes
write.csv(diff_206778, "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/preprocessed/diff_206778_DE.csv", row.names = FALSE)
# GSE288213 — normalized counts
write.csv(norm_counts_288213, "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/preprocessed/norm_counts_288213.csv")
# GSE209742 — HVG expression
write.csv(expr_209742, "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/preprocessed/expr_209742_HVG.csv")

###############################################################
# TRACK A: Differential Expression & Pathway Analysis
###############################################################

# -----------------------------
# Load required packages
# -----------------------------
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Bioconductor packages
BiocManager::install(c("DESeq2", "AnnotationDbi", "org.Mm.eg.db",
                       "clusterProfiler", "ReactomePA"))

# Load libraries
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(ReactomePA)

# -----------------------------
# 5a. Define sample metadata for DESeq2 analysis
# -----------------------------
col_data <- data.frame(
  row.names = colnames(norm_counts_288213),
  Condition = c(rep("WT", 15), rep("KO", 15))   # Adjust as needed
)

dds <- DESeqDataSetFromMatrix(countData = counts_288213,
                              colData = col_data,
                              design = ~Condition)

# Filter low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res)
res <- res[order(res$padj), ]  # sort by adjusted p-value

# Export DE results
write.csv(res, "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/results/trackA_DE_results.csv")

# -----------------------------
# 5b. Annotate DE genes (already gene symbols)
# -----------------------------
res$gene_symbol <- rownames(res)

# -----------------------------
# 5c. Pathway analysis (Reactome)
# -----------------------------
# Select significant genes
sig_genes <- rownames(res)[which(res$padj < 0.05)]

# Convert gene symbols to Entrez IDs (required for ReactomePA)
sig_genes_entrez <- mapIds(org.Mm.eg.db,
                           keys = sig_genes,
                           column = "ENTREZID",
                           keytype = "SYMBOL",
                           multiVals = "first")
sig_genes_entrez <- na.omit(sig_genes_entrez)  # remove unmapped genes

# Run Reactome enrichment
if(length(sig_genes_entrez) > 0){
  reactome_res <- enrichPathway(gene = sig_genes_entrez,
                                organism = "mouse",
                                pvalueCutoff = 0.05,
                                readable = TRUE)
  # Export pathway results
  write.csv(as.data.frame(reactome_res),
            "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/results/trackA_Reactome_results.csv")
}

# -----------------------------
# 5d. Optional: GO enrichment (Biological Process)
# -----------------------------
if(length(sig_genes_entrez) > 0){
  go_res <- enrichGO(gene          = sig_genes_entrez,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENTREZID",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     readable      = TRUE)
  
  # Export GO results
  write.csv(as.data.frame(go_res),
            "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/results/trackA_GO_BP_results.csv")
}

library(tidyverse)

# -------------------------------------------
# Paths
# -------------------------------------------
proj_dir <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA"
pre_dir  <- file.path(proj_dir, "preprocessed")
res_dir  <- file.path(proj_dir, "results")
dir.create(res_dir, showWarnings = FALSE)

# -------------------------------------------
# Load DEGs
# -------------------------------------------
deg_206778 <- read.csv(file.path(pre_dir, "diff_206778_DE.csv"))
deg_288213 <- read.csv(file.path(pre_dir, "DEGs_288213.csv"))

# -------------------------------------------
# Standardize gene column name
# -------------------------------------------
colnames(deg_206778)[1] <- "Gene"
colnames(deg_288213)[1] <- "Gene"

# -------------------------------------------
# Filter significant genes
# (keeps only genes with q_value / padj < 0.05)
# -------------------------------------------

# dataset 206778 (already has q_value)
deg_206778_sig <- deg_206778 %>% filter(!is.na(q_value) & q_value < 0.05)

# dataset 288213 (DESeq2 output uses padj)
deg_288213_sig <- deg_288213 %>% filter(!is.na(padj) & padj < 0.05)

# -------------------------------------------
# Merge on Gene
# -------------------------------------------
deg_merged <- full_join(
  deg_206778_sig,
  deg_288213_sig,
  by = "Gene",
  suffix = c("_206778", "_288213")
)

# -------------------------------------------
# Save merged DEGs
# -------------------------------------------
write.csv(deg_206778_sig, file.path(res_dir, "DEGs_206778.csv"), row.names = FALSE)
write.csv(deg_288213_sig, file.path(res_dir, "DEGs_288213_sig.csv"), row.names = FALSE)

write.csv(deg_merged, file.path(res_dir, "DEGs_merged.csv"), row.names = FALSE)

cat("Merged DEGs saved at:\n",
    file.path(res_dir, "DEGs_merged.csv"), "\n")


###############################################################
# TRACK A: Visualization & Summary of DE & Pathway Results
###############################################################

# -----------------------------
# Load required libraries
# -----------------------------
library(ggplot2)
library(pheatmap)
library(ReactomePA)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(enrichplot)

# -----------------------------
# Define results folder
# -----------------------------
results_dir <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA/results/"
if(!dir.exists(results_dir)) dir.create(results_dir)

# -----------------------------
# 1. Volcano Plot of DE genes
# -----------------------------
res$log10padj <- -log10(res$padj)
res$Significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, "Yes", "No")

volcano_plot <- ggplot(res, aes(x=log2FoldChange, y=log10padj, color=Significant)) +
  geom_point(alpha=0.5) +
  scale_color_manual(values = c("No"="grey", "Yes"="red")) +
  theme_minimal() +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  geom_hline(yintercept=-log10(0.05), col="blue", linetype="dashed") +
  geom_vline(xintercept=c(-1,1), col="blue", linetype="dashed")

ggsave(filename = paste0(results_dir, "trackA_volcano.png"),
       plot = volcano_plot, width = 7, height = 6)

# -----------------------------
# 2. Heatmap of top DE genes
# -----------------------------
top_genes <- head(rownames(res[order(res$padj), ]), 50)  # top 50 DE genes
norm_counts_top <- log2(counts(dds)[top_genes, ] + 1)

heatmap_file <- paste0(results_dir, "trackA_heatmap_top50.png")
pheatmap(norm_counts_top,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         filename = heatmap_file,
         main = "Top 50 DE Genes Heatmap")

# -----------------------------
# Reactome enrichment plots
# -----------------------------
if(exists("reactome_res") && nrow(as.data.frame(reactome_res)) > 0){
  # Barplot
  reactome_bar <- barplot(reactome_res, showCategory = 20)
  png(paste0(results_dir, "trackA_Reactome_barplot.png"), width=1000, height=800)
  print(reactome_bar)  # <- important!
  dev.off()
  
  # Dotplot
  reactome_dot <- dotplot(reactome_res, showCategory = 20)
  png(paste0(results_dir, "trackA_Reactome_dotplot.png"), width=1000, height=800)
  print(reactome_dot)  # <- important!
  dev.off()
}

# -----------------------------
# GO BP enrichment plots
# -----------------------------
if(exists("go_res") && nrow(as.data.frame(go_res)) > 0){
  # Barplot
  go_bar <- barplot(go_res, showCategory = 20)
  png(paste0(results_dir, "trackA_GO_BP_barplot.png"), width=1000, height=800)
  print(go_bar)  # <- important!
  dev.off()
  
  # Dotplot
  go_dot <- dotplot(go_res, showCategory = 20)
  png(paste0(results_dir, "trackA_GO_BP_dotplot.png"), width=1000, height=800)
  print(go_dot)  # <- important!
  dev.off()
}

###############################################################
# END OF TRACK A
###############################################################

# -----------------------------
# 5. Summary tables
# -----------------------------
# Save top 20 DE genes
write.csv(head(res[order(res$padj), ], 20),
          file = paste0(results_dir, "trackA_top20_DE_genes.csv"))

# Save top 20 Reactome pathways
if(exists("reactome_res") && nrow(as.data.frame(reactome_res)) > 0){
  write.csv(head(as.data.frame(reactome_res), 20),
            file = paste0(results_dir, "trackA_top20_Reactome.csv"))
}

# Save top 20 GO BP terms
if(exists("go_res") && nrow(as.data.frame(go_res)) > 0){
  write.csv(head(as.data.frame(go_res), 20),
            file = paste0(results_dir, "trackA_top20_GO_BP.csv"))
}
