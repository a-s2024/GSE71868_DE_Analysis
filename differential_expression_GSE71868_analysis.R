# =============================================================
# Title: Differential Expression Analysis of GSE71868 Microarray Data
# Author: a-s2024
# Description:
#   This script performs differential gene expression analysis 
#   on the GEO dataset GSE71868 using the limma package.
#   It compares gene expression in pulmonary CD11c+ cells between 
#   middle-aged and young mice to identify genes involved in age-related 
#   inflammation, immune response, and arachidonic acid metabolism.
#   Data were generated using the Illumina MouseRef-8 v2.0 expression beadchip.
# =============================================================

# ------------------ Load Required Libraries ------------------
library(Biobase)     # For working with ExpressionSet objects
library(GEOquery)    # To download and parse GEO data
library(limma)       # Linear modeling for differential expression


# ------------------ Set Up Project Paths ------------------
# Use current working directory and create a results folder
basedir <- getwd()
resfolder <- file.path(basedir, "results")
dir.create(resfolder, showWarnings = FALSE)


# ------------------ Download and Load GEO Dataset ------------------
# Download GEO dataset GSE71868 (mouse pulmonary CD11c+ cells: aging study)

gset <- getGEO("GSE71868", destdir = resfolder, GSEMatrix = TRUE)

# Select appropriate dataset if multiple found
idx <- grep("GSE71868", attr(gset, "names"))
gset <- gset[[idx]]  # Now we have the ExpressionSet object


# ------------------ Explore Dataset Structure ------------------
# Check expression matrix and sample metadata
cat("ExpressionSet Class: ", class(gset), "\n")
cat("Number of samples:", ncol(exprs(gset)), "\n")
cat("Number of genes:", nrow(exprs(gset)), "\n")

# View sample and gene-level annotations
head(pData(gset))
head(fData(gset))


# ------------------ Prepare Metadata for Analysis ------------------
# Ensure feature variable labels are syntactically valid
fvarLabels(gset) <- make.names(fvarLabels(gset))

# Define sample groups: 4 Middle-aged vs. 4 Young samples
group_labels <- c("Middle_aged", "Middle_aged", "Middle_aged", "Middle_aged",
                  "Young", "Young", "Young", "Young")
group_factor <- as.factor(group_labels)

# Assign group labels to dataset
gset$description <- group_factor

# Create design matrix for limma analysis
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(group_factor)


# ------------------ Differential Expression Analysis ------------------
# Fit linear model to expression data
fit <- lmFit(gset, design)

# Define comparison of interest: Middle-aged vs. Young
cont.matrix <- makeContrasts(Middle_aged - Young, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)  # Empirical Bayes moderation

# Extract ranked table of differentially expressed genes
tT <- topTable(fit2, adjust = "fdr", sort.by = "B", number = nrow(fit2))


# ------------------ Filter and View Significant Genes ------------------
# Identify genes with adjusted p-value < 0.01
significant_genes <- tT$ID[tT$adj.P.Val < 0.01]
cat("Number of significantly differentially expressed genes:", length(significant_genes), "\n")
print(significant_genes)

# Save full results to CSV
write.csv(tT, file.path(resfolder, "differential_expression_results.csv"))


# ------------------ Volcano Plot Visualization ------------------
# Visualize fold changes vs. significance
plot(tT$logFC, 1 - tT$adj.P.Val,
     xlim = c(-6, 6),
     main = "Volcano Plot of Differentially Expressed Genes",
     xlab = "Log2 Fold Change",
     ylab = "1 - Adjusted P-value")
abline(h = 0.95, col = "red")  # Highlight threshold line


# ------------------ Hierarchical Clustering ------------------
# Cluster top significant genes for visual inspection
selected_genes <- head(tT$ID)  # Take top 6 genes
exprs_subset <- exprs(gset)[selected_genes, ]
exprs_log <- log2(exprs_subset)  # Log-transform expression values

# Compute distance and perform clustering
d <- dist(exprs_log)
hc <- hclust(d)

# Plot dendrogram
plot(hc, main = "Dendrogram of Top Differentially Expressed Genes")
