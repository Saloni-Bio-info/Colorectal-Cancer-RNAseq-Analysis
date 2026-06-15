# Colorectal Cancer RNA-seq Analysis (GSE104836)

## Overview

This project investigates transcriptomic alterations in colorectal cancer using RNA-seq data from the GEO dataset GSE104836. The analysis integrates differential expression analysis, functional enrichment, protein–protein interaction (PPI) network construction, and machine learning-based biomarker selection to identify potential molecular signatures associated with colorectal cancer.

## Workflow

1. Differential Gene Expression (DEG) Analysis
2. Principal Component Analysis (PCA)
3. Volcano Plot and Heatmap Visualization
4. Gene Ontology (GO) Enrichment Analysis
5. KEGG Pathway Enrichment Analysis
6. Protein–Protein Interaction (PPI) Network Construction
7. Hub Gene Identification using Cytoscape
8. LASSO-Based Biomarker Selection
9. ROC Curve Evaluation

---

## Dataset

* **GEO Accession:** GSE104836
* **Data Source:** NCBI Gene Expression Omnibus (GEO)
* **Organism:** Homo sapiens
* **Sample Types:** Colorectal cancer and normal tissue samples

---

## Methods

### 1. Differential Expression Analysis

* RNA-seq count data processing
* DESeq2 normalization and statistical testing
* Identification of significantly differentially expressed genes (DEGs)

### 2. Visualization

* Principal Component Analysis (PCA)
* Volcano Plot
* Heatmap of significant DEGs

### 3. Functional Enrichment Analysis

* Gene Ontology (GO) enrichment
* KEGG pathway enrichment

### 4. PPI Network Analysis

* STRING database interaction network
* Cytoscape visualization
* Network topology analysis
* Hub gene identification using Degree Centrality

### 5. Machine Learning-Based Biomarker Selection

* LASSO logistic regression using hub genes
* Feature selection for candidate biomarker discovery
* ROC curve analysis for model evaluation

---

## Repository Structure

```text
Data/            - Input expression data
DEGs/            - Differentially expressed genes
Enrichment/      - GO and KEGG enrichment results
Figures/         - PCA, Heatmap, Volcano Plot
ML Biomarker/    - LASSO biomarker selection results, ROC Curve
PPI/             - STRING and Cytoscape network analysis
updatedAnalysis.R - Main analysis script
```

---

## Key Outcomes

### Differential Expression Analysis

* Identification of significantly dysregulated genes in colorectal cancer samples.

### Functional Enrichment

* Detection of enriched biological processes and signaling pathways associated with colorectal cancer progression.

### PPI Network Analysis

* Construction of a colorectal cancer interaction network.
* Identification of hub genes with high network connectivity.

### Machine Learning Biomarker Selection

LASSO logistic regression identified the following candidate biomarkers:

* MMP1
* COL1A1
* MMP3
* MMP7
* COL11A1
* ACAN
* INHBA

### Model Performance

The biomarker panel was evaluated using Receiver Operating Characteristic (ROC) analysis and achieved:

* **AUC = 1.00**

indicating perfect separation of colorectal cancer and normal samples within the GSE104836 dataset.

---

## Limitations

* The dataset contains only 20 samples.
* Model evaluation was performed on the same dataset used for training.
* Independent validation using additional colorectal cancer datasets is required to assess generalizability and clinical utility.

---

## Future Directions

* External validation using independent GEO datasets
* Comparison with Random Forest and XGBoost models
* Survival analysis of identified biomarkers
* Integration of multi-omics datasets for biomarker discovery

---

## Tools and Packages

* R
* DESeq2
* clusterProfiler
* STRING
* Cytoscape
* glmnet
* pROC
