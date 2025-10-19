#####################################################################
# HasmikChilingaryan_ClassIII_DEG.R
# Module III — Class 3C: Microarray Data Analysis (GSE1977)
# Workflow: Probe IDs to gene mapping, Differential Gene Expression
# Analysis, Data Visualization
# Working Data from:
# Module_III_GSE1977_DEG\raw_data\GSE1977_series_matrix.txt
#####################################################################

# housekeeping
rm(list = ls()); gc()
options(stringsAsFactors = FALSE)
getwd()
setwd("\\\\adf\\storage\\h\\c\\hxc411\\AI_BiotechBioinf_Internship2025\\Module_III_GSE1977_DEG")
list.files()
list.files("raw-data")

# insatll and load packages
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "AnnotationDbi", "hgu95av2.db"), ask = FALSE)
install.packages(c("dplyr", "tibble", "ggplot2", "pheatmap"), dependencies = TRUE)

# load libraries
library(GEOquery)
library(Biobase)
library(limma)
library(AnnotationDbi)
library(hgu95av2.db)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)

# load GEO series matrix
series_file <- "raw_data/GSE1977_series_matrix.txt"
library(GEOquery)
gse <- getGEO(filename = series_file)
class(gse) # should be list
length(gse) # number of expressionsets
eset <- gse

# extract expression matrix and phenotype
exprs_matrix <- exprs(eset)
pheno <- pData(eset)

cat("Expression matrix dimensions (probes x samples):", dim(exprs_matrix), "\n")
cat("Phenotype dimensions:", dim(pheno), "\n")

# Probe-to-Gene Mapping
probe_ids <- rownames(exprs_matrix)

# map probe IDs to gene symbols
gene_symbols <- mapIds(
  hgu95av2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

# merge gene symbols with expression matrix
exprs_df <- as.data.frame(exprs_matrix)
exprs_df$SYMBOL <- gene_symbols[rownames(exprs_df)]
exprs_df <- exprs_df[!is.na(exprs_df$SYMBOL), ]  # remove probes with no symbol

# collapse multiple probes per gene by averaging
averaged_exprs <- limma::avereps(as.matrix(exprs_df[,1:(ncol(exprs_df)-1)]), ID = exprs_df$SYMBOL)

cat("Expression matrix after mapping and averaging duplicates: ", dim(averaged_exprs), "\n")

# Define sample groups (IR, UV, Mock)
# using "title" column in phenotype for treatment
pheno$group <- sapply(pheno$title, function(x){
  if(grepl("IR", x)) return("IR")
  else if(grepl("UV", x)) return("UV")
  else return("Mock")
})
groups <- factor(pheno$group)
levels(groups)

# Differential Gene Expression Analysis using limma
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups)

# example contrast: IR vs Mock
contrast_matrix <- makeContrasts(IR_vs_Mock = IR - Mock,
                                 UV_vs_Mock = UV - Mock,
                                 levels = design)

fit <- lmFit(averaged_exprs, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# extract DEGs for IR vs Mock
deg_results <- topTable(fit2, coef = "IR_vs_Mock", number = Inf, adjust.method = "BH")

# classify DEGs
deg_results$threshold <- as.factor(ifelse(
  deg_results$adj.P.Val < 0.05 & deg_results$logFC > 1, "Upregulated",
  ifelse(deg_results$adj.P.Val < 0.05 & deg_results$logFC < -1, "Downregulated", "No")
))

# split by regulation
upregulated <- subset(deg_results, threshold == "Upregulated")
downregulated <- subset(deg_results, threshold == "Downregulated")
deg_updown <- rbind(upregulated, downregulated)

# save DEG tables
write.csv(deg_results, file = "results/DEGs_Results.csv")
write.csv(upregulated, file = "results/Upregulated_DEGs.csv")
write.csv(downregulated, file = "results/Downregulated_DEGs.csv")
write.csv(deg_updown, file = "results/Updown_DEGs.csv")

# Data Visualization

# volcano plot
png("plots/volcano_plot.png", width = 2000, height = 1500, res = 300)
ggplot(deg_results, aes(x = logFC, y = -log10(adj.P.Val), color = threshold)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "No" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of Differentially Expressed Genes (IR vs Mock)",
       x = "log2 Fold Change",
       y = "-log10(Adjusted P-value)",
       color = "Regulation")
dev.off()

# heatmap of top 25 DEGs
top_genes <- head(rownames(deg_updown[order(deg_updown$adj.P.Val), ]), 25)
heatmap_data <- averaged_exprs[top_genes, ]

# column names formatted by group
group_char <- as.character(groups)
heatmap_names <- ave(group_char, group_char, FUN = function(x) paste0(x, "_", seq_along(x)))
colnames(heatmap_data) <- heatmap_names

png("plots/heatmap_top25_DEGs.png", width = 2000, height = 1500, res = 300)
pheatmap(
  heatmap_data,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 6,
  fontsize_col = 8,
  main = "Top 25 Differentially Expressed Genes (IR vs Mock)"
)
dev.off()

# Assignment Summary (4–5 lines)

cat("\nAssignment Summary:\n")
cat("1. Multiple probes mapping to the same gene were averaged using limma::avereps().\n")
cat("2. Contrast performed: IR vs Mock.\n")
cat(paste0("3. Upregulated genes: ", nrow(upregulated), "\n"))
cat(paste0("4. Downregulated genes: ", nrow(downregulated), "\n"))
cat("5. DEG tables and plots saved in 'results/' and 'plots/' folders.\n")