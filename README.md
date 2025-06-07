# Transcriptome-Profiling-of-Proliferating-Senescent-and-Mitochondrially-Depleted-Cells-in-R
A comprehensive R-based workflow for differential expression, visualization, and functional enrichment analysis across proliferating, senescent, and mitochondrially-depleted cell states, developed as coursework for MSc Bioinformatics 2024/25, University of Glasgow.

# Overview
This repository provides a comprehensive workflow for analyzing bulk RNA-seq data from human fibroblasts (IMR90) under proliferating, senescent, and mitochondria-depleted senescent conditions.

The `omic_functions.R` script contains a suite of custom functions that handle data loading, preprocessing, differential expression analysis, visualization (including volcano and MA plots, PCA, heatmaps), and gene enrichment analyses such as Gene Ontology and Gene Set Enrichment Analysis. These functions enable streamlined, reproducible analysis tailored to RNA-seq datasets.

The `main_script.R` script applies these functions in a structured pipeline. It loads expression data, differential expression results, and annotations, integrates and cleans the data, performs statistical analyses, generates multiple visualizations, identifies significant gene signatures, and runs functional enrichment analyses. This script coordinates the entire analysis workflow, producing interpretable results and figures to characterize transcriptional changes across the different cell states.

Together, the scripts support detailed exploration and interpretation of transcriptional differences relevant to cellular senescence biology.


# omic\_functions.R — RNA-seq Analysis Functions

## Features

* **Data Import and Preprocessing**
  Functions for reading CSV expression data, merging expression matrices with differential expression (DE) results, and annotation data, generating master tables for downstream analysis.

* **Differential Expression Analysis Helpers**
  Functions to create filtered tables of significantly differentially expressed genes (adjusted p-value < 0.001 and |log2FC| > 2), merge DE tables from multiple comparisons, and generate tables keyed by gene symbols.

* **Plotting and Visualization**

  * Volcano and MA plots for DE visualization
  * PCA and density plots for quality and sample comparison
  * Boxplots, violin plots, jitter plots for gene expression visualization, including combined violin-jitter plots
  * Fold-vs-Fold comparison plots between conditions
  * Clustered heatmaps for expression patterns of gene sets
  * Venn diagrams and overlap significance tests (hypergeometric tests)

* **Enrichment and Network Analysis**

  * Over-Representation Analysis (ORA) and Gene Set Enrichment Analysis (GSEA) using `clusterProfiler`
  * Functional plots such as barplots, dotplots, cnetplots, goplots, and ridgeplots
  * STRING database integration for protein-protein interaction network visualization
  * Custom metagene boxplots summarizing expression across gene signatures

* **Utility and Theming**

  * Custom ggplot2 themes for consistent figure styling
  * Plot saving functions with customizable file paths and sizes

---

## Usage

* Source this script in your main analysis script using:

  ```r
  source("omic_functions.R")
  ```
* The functions can be used to perform data wrangling, visualization, and biological interpretation of RNA-seq results.
* Many functions require inputs such as expression matrices, DE result tables, annotation data, and parameters specifying sample groups and colors.
* The script contains functions to facilitate a workflow focusing on senescence-associated transcriptomic changes, immune and cell cycle gene expression, and mitochondrial depletion effects.


## Dependencies

This script requires the following R packages:

* ggplot2
* ggrepel
* eulerr
* reshape2
* amap
* clusterProfiler
* org.Hs.eg.db
* STRINGdb

Make sure these are installed before running the functions (or uncomment the lines in the script to do so):

```r
install.packages(c("ggplot2", "ggrepel", "eulerr", "reshape2", "amap"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "STRINGdb"))
```

# main_script.R — DGE Analysis using custom functions

## Workflow Summary

### 1. Setup & Data Loading

* Load custom functions from `omic_functions.R`.
* Read expression matrix (EM), differential expression (DE) results for three pairwise comparisons, gene annotations, and sample metadata.
* Use `csv_reader()` custom function to load CSVs efficiently.

### 2. Master Table Creation

* For each pairwise DE comparison, create a "master" table by combining expression, DE results, and annotations (`create_master_table()`).
* Merge the master tables for all comparisons to a single master table for integrated analysis.
* Clean and reorder columns for clarity.
* Order genes by minimum adjusted p-value across all comparisons.

### 3. Expression Data Processing

* Create a gene expression table indexed by gene symbols (`create_symbols_table()`).
* Scale expression data to normalize for PCA and clustering.

### 4. Subset Significant Genes

* Define strict thresholds (|log2FC| > 2 and adjusted p < 0.001) to identify:

  * Genes significant in **all three comparisons**
  * Genes significant in **any comparison** (commented out here, but easy to enable)
* Extract gene symbol vectors for downstream analyses.

### 5. Visualization

* Generate **volcano plots** and **MA plots** per comparison to visualize DEGs.
* Plot **expression density** to check data distribution.
* Produce **PCA plots** using scaled expression to visualize sample clustering by condition.
* Create **boxplots** for top significant genes per comparison.
* Generate **clustered heatmaps** of significant genes.
* Plot **rugplots** to annotate sample groups on heatmaps.
* (Optional) Create Venn diagrams to visualize gene overlaps.

### 6. Gene Overlap & Correlation

* Calculate hypergeometric tests to evaluate overlap significance between DEG sets.
* Perform correlation tests for fold changes between pairwise comparisons.
* Generate fold-change vs fold-change (FVF) plots (optional, commented out).

### 7. Upregulated and Downregulated Gene Subsets

* Identify genes consistently upregulated or downregulated across all comparisons for focused analysis.

### 8. Gene ID Conversion

* Convert gene symbols to ENTREZ IDs for enrichment analyses using `bitr()` from `clusterProfiler`.
* Conversion supports further ORA (Overrepresentation Analysis) and GSEA.

### 9. Functional Enrichment Analyses

* Perform **GO enrichment (ORA)** for biological processes (BP) on:

  * All significant genes
  * Upregulated genes
  * Downregulated genes
* Perform **GSEA** for each pairwise comparison (optional).
* Visualize enrichment results with barplots and clustered heatmaps.

### 10. Signature Analysis

* Define custom gene signatures based on combined fold change and significance criteria (e.g., upregulated in one comparison but downregulated in another).
* Convert signature genes to ENTREZ IDs.
* Perform functional enrichment (ORA) on signature gene sets.
* Generate signature-specific boxplots and heatmaps.

## Tips for Customization

* **Comparisons:** Change which DE results to analyze by modifying `create_master_table()` inputs or filtering criteria.
* **Signatures:** Define your own signatures with logical expressions on log2FC and significance columns in the master table.
* **Thresholds:** Adjust |log2FC| and p-adj cutoffs as needed.
* **Plots:** Modify `plot_volcano()`, `plot_ma()`, `plot_pca()`, and other visualization functions parameters for titles, file names, or aesthetics.
* **Enrichment:** Expand enrichment analysis to other ontologies (MF, CC) or gene sets by modifying `do_ora()` and `do_gsea()` calls.
* **Save paths:** Ensure all generated plots and tables are saved with clear file names to prevent overwriting.

