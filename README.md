# Differential Gene Expression Analysis: GSE71868

This project performs a differential expression analysis on the GSE71868 dataset, which compares gene expression between pulmonary CD11c+ cells from young and middle-aged mice.

## Objective

Identify genes that are significantly differentially expressed with aging in pulmonary CD11c+ cells using the limma package in R.

## Dataset

- Source: [GEO - GSE71868](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71868)
- Organism: Mus musculus (mouse)
- Platform: GPL6885 — Illumina MouseRef-8 v2.0 expression beadchip
- Samples: 4 young mice, 4 middle-aged mice

## Tools Used

- R
- GEOquery
- limma
- Biobase

## Output

- `differential_expression_results.csv`: Ranked list of genes with logFC and adjusted p-values
- Volcano plot and dendrogram

## Structure
- `gene_expression_analysis.R` *Main R script performing the analysis*
- `differential_expression_results.csv` *Output CSV with significant genes and stats*
- `volcano_plot.png` *Volcano plot visualization of DE genes*
- `dendrogram.png` *Dendrogram of hierarchical clustering*
- `README.md` *This file*
