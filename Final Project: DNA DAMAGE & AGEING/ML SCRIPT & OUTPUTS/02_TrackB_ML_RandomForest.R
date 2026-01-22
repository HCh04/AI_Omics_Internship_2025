###############################################################
# TRACK B — MACHINE LEARNING FRAMEWORK
# Final Project: DNA Damage & Ageing
###############################################################

# -----------------------------
# Load required packages
# -----------------------------
library(tidyverse)
library(caret)
library(randomForest)
library(ggplot2)

# -----------------------------
# Define paths
# -----------------------------
proj_dir  <- "U:/AI_BiotechBioinf_Internship2025/Final_Project_DDA"
pre_dir   <- file.path(proj_dir, "preprocessed")
raw_dir   <- file.path(proj_dir, "data_raw")
res_dir   <- file.path(proj_dir, "results")
dir.create(res_dir, showWarnings = FALSE)

# -----------------------------
# Load expression matrix
# -----------------------------
expr_file <- file.path(pre_dir, "norm_counts_288213.csv")
expr <- read.csv(expr_file, row.names = 1)
expr <- as.data.frame(expr)

# Transpose so samples are rows
expr_t <- as.data.frame(t(expr))

# Make column names valid
colnames(expr_t) <- make.names(colnames(expr_t))

# -----------------------------
# Load cleaned metadata
# -----------------------------
meta_file <- file.path(raw_dir, "GSE288213_metadata_clean.csv")
meta <- read.csv(meta_file, stringsAsFactors = FALSE)

# -----------------------------
# Map metadata sample names to expression rownames
# -----------------------------
mapping <- data.frame(
  expr_sample = c("Csa_1m15d_01","Csa_1m15d_02","Csa_1m15d_03","Csa_1m15d_04","Csa_1m15d_05",
                  "Csa_12m_01","Csa_12m_02","Csa_12m_03","Csa_12m_04","Csa_12m_05",
                  "Csa_24m_01","Csa_24m_02","Csa_24m_03","Csa_24m_04","Csa_24m_05",
                  "WT_1m15d_01","WT_1m15d_02","WT_1m15d_03","WT_1m15d_04","WT_1m15d_05",
                  "WT_12m_01","WT_12m_02","WT_12m_03","WT_12m_04","WT_12m_05",
                  "WT_24m_01","WT_24m_02","WT_24m_03","WT_24m_04","WT_24m_05"),
  meta_sample = c("SRR32151315","SRR32151325","SRR32151326","SRR32151327","SRR32151328",
                  "SRR32151329","SRR32151330","SRR32151331","SRR32151335","SRR32151336",
                  "SRR32151337","SRR32151338","SRR32151339","SRR32151343","SRR32151316",
                  "SRR32151317","SRR32151318","SRR32151319","SRR32151320","SRR32151321",
                  "SRR32151322","SRR32151323","SRR32151324","SRR32151332","SRR32151333",
                  "SRR32151334","SRR32151340","SRR32151341","SRR32151342","SRR32151344")
)
meta$sample <- mapping$expr_sample[match(meta$sample, mapping$meta_sample)]

# -----------------------------
# Align samples
# -----------------------------
common_samples <- intersect(rownames(expr_t), meta$sample)
expr_t <- expr_t[common_samples, ]
meta <- meta[meta$sample %in% common_samples, ]
meta <- meta[match(rownames(expr_t), meta$sample), ]
stopifnot(all(rownames(expr_t) == meta$sample))
cat("Sample names aligned. Number of samples:", nrow(expr_t), "\n")

# -----------------------------
# Attach labels
# -----------------------------
expr_t$genotype <- as.factor(meta$genotype)
expr_t$age_days <- meta$age_days

# -----------------------------
# Make column names valid again
colnames(expr_t) <- make.names(colnames(expr_t))

# -----------------------------
# Select top variable genes
num_features <- 500
gene_vars <- apply(expr_t[, sapply(expr_t, is.numeric)], 2, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:num_features]

ml_data <- expr_t[, c(top_genes, "genotype")]

# -----------------------------
# Standardize numeric features
num_cols <- sapply(ml_data, is.numeric)
ml_data[, num_cols] <- scale(ml_data[, num_cols])

# -----------------------------
# Train/Test split
set.seed(123)
train_index <- createDataPartition(ml_data$genotype, p = 0.7, list = FALSE)
train <- ml_data[train_index, ]
test  <- ml_data[-train_index, ]

# -----------------------------
# Train Random Forest
rf_model <- randomForest(
  genotype ~ .,
  data = train,
  ntree = 500,
  importance = TRUE
)

# -----------------------------
# Save model
saveRDS(rf_model, file = file.path(res_dir, "rf_model.rds"))

# -----------------------------
# Model Evaluation
pred <- predict(rf_model, test)
cm <- confusionMatrix(pred, test$genotype)

# Save confusion matrix
write.csv(as.data.frame(cm$table), file = file.path(res_dir, "confusion_matrix.csv"))
cat("Confusion matrix saved.\n")
print(cm)

# -----------------------------
# Variable Importance
varImpPlot(rf_model, n.var = 20)
imp <- importance(rf_model)
imp_df <- data.frame(Gene = rownames(imp), MeanDecreaseGini = imp[, "MeanDecreaseGini"])
imp_df <- imp_df %>% arrange(desc(MeanDecreaseGini))

# Save variable importance table
write.csv(imp_df, file = file.path(res_dir, "variable_importance.csv"))
cat("Variable importance table saved.\n")

# Save top 20 variable importance plot
p <- ggplot(imp_df[1:20, ], aes(x = reorder(Gene, MeanDecreaseGini), y = MeanDecreaseGini)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 20 Important Genes", x = "Gene", y = "Mean Decrease Gini")
ggsave(filename = file.path(res_dir, "top20_variable_importance.png"), plot = p, width = 8, height = 6)
cat("Top 20 variable importance plot saved.\n")
