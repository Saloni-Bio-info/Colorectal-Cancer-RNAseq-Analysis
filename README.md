Colorectal Cancer RNA-seq Analysis

Overview

This project investigates transcriptomic alterations in colorectal cancer using RNA-seq data from the GEO database. The workflow includes differential expression analysis, functional enrichment analysis, protein-protein interaction (PPI) network construction, and hub gene identification.

Dataset

* GEO Accession: GSE104836
* Organism: Homo sapiens
* Disease: Colorectal Cancer
* Samples: 20 (10 Tumor vs 10 Normal)

Methods

Differential Expression Analysis

* DESeq2
* Count normalization
* Statistical testing
* Multiple testing correction

Functional Enrichment Analysis

* Gene Ontology (GO)
* KEGG Pathway Analysis
* clusterProfiler

Network Biology

* STRING Database
* Cytoscape
* Degree Centrality Analysis

Results

Differentially Expressed Genes

* 4,734 significant DEGs identified

Functional Enrichment Highlights

- Extracellular matrix organization
- Extracellular structure organization
- Muscle system process
- Muscle contraction
- Cell cycle pathway
- Cytokine-cytokine receptor interaction

Hub Genes Identified

* MMP1
* COL1A1
* MMP3
* MMP7
* BGN
* COL11A1
* ACAN
* ADAM12
* ADAMTS2
* COL7A1

Tools Used

* R
* RStudio
* DESeq2
* clusterProfiler
* STRING
* Cytoscape

Repository Structure

* Data/
* DEGs/
* Enrichment/
* Figures/
* PPI/
* ML_Biomarker/

Planned Future Work
- Machine Learning-based Biomarker Discovery

Author

Saloni 
M.Sc. Bioinformatics
