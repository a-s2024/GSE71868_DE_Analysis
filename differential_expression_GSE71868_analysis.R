# =============================================================
# Title: Differential Expression Analysis of GSE71868 Microarray Data
# Author: a-s2024
# Description:
#   This script performs differential gene expression analysis
#   on GEO dataset GSE71868 (human brain aging study).
# =============================================================

# ------------------ Load Required Libraries ------------------
library(Biobase)     # For ExpressionSet objects
library(GEOquery)    # Download GEO data
library(limma)       # Differential expression analysis

# ------------------ Set Up Project Paths ------------------
basedir <- getwd()
resfolder <- file.path(basedir, "results")
dir.create(resfolder, showWarnings = FALSE)

# ------------------ Download and Load GEO Dataset ------------------
gset <- getGEO("GSE71868", destdir = resfolder, GSEMatrix = TRUE)
idx <- grep("GSE71868", attr(gset, "names"))
gset <- gset[[idx]]  # Extract ExpressionSet object

# ------------------ Explore Dataset ------------------
cat("ExpressionSet Class:", class(gset), "\n")
cat("Number of samples:", ncol(exprs(gset)), "\n")
cat("Number of genes:", nrow(exprs(gset)), "\n")

# ------------------ Preprocessing: Check & Normalize ------------------
cat("Summary of expression values before any transform:\n")
print(summary(exprs(gset)[,1]))

# Optional histogram to visualize distribution
hist(exprs(gset), breaks = 100, main = "Histogram of Expression Values", xlab = "Expression")

# Check if data appears log2 transformed or not
if (min(exprs(gset)) > 0 & max(exprs(gset)) < 20) {
  cat("Data appears already log2 transformed; skipping log2 transform.\n")
} else {
  cat("Data NOT log2 transformed; applying log2 transform.\n")
  exprs(gset) <- log2(exprs(gset) + 1)  # Add 1 to avoid log(0)
}

# Normalize expression matrix (quantile normalization)
expr_matrix <- exprs(gset)
expr_norm <- normalizeBetweenArrays(expr_matrix, method = "quantile")
exprs(gset) <- expr_norm

cat("Summary of expression values after normalization:\n")
print(summary(exprs(gset)[,1]))

# ------------------ Prepare Metadata ------------------
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Define sample groups (adjust as per your samples)
group_labels <- c("Middle_aged", "Middle_aged", "Middle_aged", "Middle_aged",
                  "Young", "Young", "Young", "Young")
group_factor <- as.factor(group_labels)
gset$description <- group_factor

# Create design matrix
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(group_factor)

# ------------------ Differential Expression Analysis ------------------
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(Middle_aged - Young, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)  # Empirical Bayes moderation

# Extract ranked table of genes, sorted by B statistic
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = nrow(fit2))
summary(tT$logFC)

# ------------------ Significant Genes ------------------
significant_genes <- tT$ID[tT$adj.P.Val < 0.01]
cat("Number of significantly differentially expressed genes:", length(significant_genes), "\n")
print(significant_genes)

# Save results
write.csv(tT, file.path(resfolder, "differential_expression_results.csv"))

# ------------------ Volcano Plot ------------------
plot(tT$logFC, 1 - tT$adj.P.Val,
     xlim = c(-6, 6),
     main = "Volcano Plot of Differentially Expressed Genes",
     xlab = "Log2 Fold Change",
     ylab = "1 - Adjusted P-value")
abline(h = 0.95, col = "red")  # Threshold line

# ------------------ Hierarchical Clustering ------------------
# Select top 10 significant genes
selected_genes <- head(tT$ID, 10)

# Subset expression data for those genes
exprs_subset <- exprs(gset)[selected_genes, ]

# Make sure rownames are gene names for labeling
rownames(exprs_subset) <- selected_genes

# Compute Euclidean distance between genes
d <- dist(exprs_subset)

# Perform hierarchical clustering using complete linkage
hc <- hclust(d, method = "complete")

# Plot dendrogram with gene names as labels
plot(hc, 
     main = "Dendrogram of Top 10 Differentially Expressed Genes",
     xlab = "Genes",
     ylab = "Distance (Euclidean)",
     sub = "",
     hang = -1)

head(tT$ID, 10)
