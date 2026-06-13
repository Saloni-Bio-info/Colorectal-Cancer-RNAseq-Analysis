Colorectal Cancer RNA-seq Analysis (GSE104836)

Overview

This project investigates transcriptomic alterations in colorectal cancer using RNA-seq data from GEO dataset GSE104836.

The workflow includes:

* Differential Gene Expression (DEG) analysis
* Functional enrichment analysis (GO and KEGG)
* Protein–Protein Interaction (PPI) network construction
* Hub gene identification using Cytoscape

Dataset

* GEO Accession: GSE104836
* Data Source: NCBI GEO
* Organism: Homo sapiens

Methods

1. Differential Expression Analysis

* R
* DESeq2
* Filtering of significant DEGs

2. Visualization

* PCA
* Volcano Plot
* Heatmap

3. Functional Enrichment

* Gene Ontology (GO)
* KEGG Pathway Analysis

4. PPI Network Analysis

* STRING database
* Cytoscape
* Network topology analysis
* Hub gene identification using Degree Centrality

Repository Structure

Data/ - Input expression data

DEGs/ - Differentially expressed genes

Enrichment/ - GO and KEGG enrichment results

Figures/ - PCA, Heatmap, Volcano Plot

PPI/ - STRING and Cytoscape network analysis

Analysis.R - Main analysis script

Key Outcomes

* Identification of significantly dysregulated genes
* Functional enrichment of biological pathways
* Construction of a colorectal cancer PPI network
* Identification of candidate hub genes for further biomarker investigation

Future Work

* Machine learning-based biomarker discovery
* Feature selection and model development
* Multi-omics integration
