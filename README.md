# Single-cell RNA-seq Analysis of Human Pancreatic Islets

## Overview

This repository contains a reproducible pipeline for analyzing human pancreatic islet single-cell RNA-seq data. The workflow includes:

* Data download (GEO: GSExxxxx)
* Quality control and filtering
* Normalization (SCTransform)
* Dimensionality reduction (PCA, UMAP)
* Clustering and visualization
* Cell type annotation
* Endocrine subclustering
* Beta cell state identification
* Differential expression analysis
* Gene Ontology enrichment

---

## Project Structure

```
islet_scRNA_pipeline/
│
├── data/
│   ├── raw/
│   ├── processed/
│
├── scripts/
│   ├── 01_download_data.R
│   ├── 02_load_and_merge.R
│   ├── 03_QC_and_clustering.R
│   ├── 04_annotation_and_states.R
│   ├── 05_DE_analysis.R
│   ├── run_pipeline.R
│
├── plots/
├── results/
├── README.md
```

---

## Environment Setup (Conda)

### Step 1: Create conda environment

```
conda create -n islet_scRNA_env r-base=4.2 -y
```

---

### Step 2: Activate environment

```
conda activate islet_scRNA_env
```

---

### Step 3: Install required system libraries

```
conda install -c conda-forge r-essentials r-devtools -y
```

---

## Install R packages

Open R inside the conda environment:

```
R
```

Then run:

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork"))
install.packages("sctransform")

install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
```

---

## Running the Pipeline

### Run step-by-step (recommended first run)

```r
source("scripts/03_QC_and_clustering.R")
source("scripts/04_annotation_and_states.R")
source("scripts/05_DE_analysis.R")
```

---

### Run full pipeline (one command)

```r
source("scripts/run_pipeline.R")
```

---

### Optional: Run from terminal

```
Rscript scripts/run_pipeline.R
```

---

## Output

### Processed data:

```
data/processed/
```

### Plots:

```
plots/
results/
```

### Key outputs:

* UMAP plots
* Cluster markers
* Beta cell state markers
* Differential expression tables
* GO enrichment results

---

## Notes

* Ensure input file exists:

```
data/processed/islet_merged.rds
```

* Pipeline assumes human gene symbols (required for GO analysis)

---

## Reproducibility

For reproducibility, the pipeline:

* Uses fixed QC thresholds
* Saves intermediate objects
* Separates preprocessing and biological analysis

---

## Citation / Use

If you use this pipeline, please cite relevant tools:

* Seurat
* clusterProfiler

---

## Author

Nishant Kumar
Computational Biology | Single-cell Analysis
