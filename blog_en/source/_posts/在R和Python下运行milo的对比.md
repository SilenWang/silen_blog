---
title: Comparing Milo Implementations in R and Python
categories: Bioinformatics
date: 2025-10-28 15:57:51
tags: [Single-cell analysis, Differential abundance analysis, Milo, Bioinformatics]
---

Milo is a differential abundance analysis method for single-cell RNA sequencing data that can detect compositional changes in cell neighborhoods across different conditions.
<!-- more -->

## Introduction to the Milo Algorithm

[Milo](https://www.nature.com/articles/s41587-021-01033-z) is a differential abundance analysis method specifically designed for single-cell RNA sequencing data. Its core idea is to detect changes in cell population composition across different experimental conditions by constructing cell neighborhoods.

My understanding is that the purpose of this algorithm is to identify cells with population proportion differences between two conditions without pre-defining them. Once detection results are obtained, these differentially abundant cells can be extracted for further analysis.

## Milo Implementation in R Environment

In R language, the miloR package is used for analysis. Here's the corresponding implementation:

```r
# Load required libraries
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(tidyr)

# Read Seurat object and convert to SingleCellExperiment
seurat_obj <- readRDS("seurat.rds")
sce <- as.SingleCellExperiment(seurat_obj)

# Set random seed for reproducibility - mainly affects the kNN component
set.seed(4466)

# Create Milo object
milo_obj <- Milo(sce)

# Build graph structure - note: k is for kNN, d is for dimensions/components
milo_obj <- buildGraph(
    milo_obj,
    k = 30, d = 20,
    transposed = TRUE, 
    reduced.dim = "UMAP"
)

# Define neighborhoods - note: k, d, and reduced_dims should match previous step
# refinement_scheme enables a faster algorithm
# Provides significant computational acceleration for large datasets
milo_obj <- makeNhoods(
    x = milo_obj,
    prop = 0.1, k = 30, d = 20,
    refined = TRUE,
    reduced_dims = "UMAP",
    refinement_scheme = "graph"
)

# Count cells in neighborhoods
milo_obj <- countCells(
    milo_obj,
    meta.data = data.frame(colData(milo_obj)),
    sample = "sample_id"
)

# Prepare experimental design - create a data frame with sample IDs and group tags
traj_design <- data.frame(colData(milo_obj))[,c("sample_id", "group")] %>%
    distinct() %>%
    mutate(
        sample_id = as.character(sample_id),
        group = as.character(group)
    )
rownames(traj_design) <- traj_design$sample_id
traj_design <- traj_design[colnames(nhoodCounts(milo_obj)), , drop=FALSE]

# Differential abundance testing - fdr.weighting = "graph-overlap"
# Works with parameters from makeNhoods to improve computation speed
# Note: Comparison groups can be specified using formula-like syntax
# design = ~ 0 + group means using group variable for grouping, ~ 0 + is fixed syntax
# model.contrasts = c('groupCase - groupControl') means:
# Use rows where group equals 'Case' as target, and 'Control' as reference
da_results <- testNhoods(
    milo_obj,
    design = ~ 0 + group,
    model.contrasts = c('groupCase - groupControl'),
    design.df = traj_design,
    fdr.weighting = "graph-overlap"
)

# View significant results
da_results %>%
  arrange(SpatialFDR) %>%
  filter(SpatialFDR < 0.1) %>%
  tail(3)

# Build neighborhood graph
milo_obj <- buildNhoodGraph(milo_obj)

# All computation steps completed
# Subsequent visualization can use plotUMAP or other functions
```

## Milo Implementation in Python Environment

Python has multiple Milo implementations. Besides [Mipopy](https://github.com/emdann/milopy) mentioned in the official [miloR project](https://github.com/MarioniLab/miloR), the [pertpy](https://github.com/scverse/pertpy) library also implements the Milo algorithm. Here's an example using pertpy:

```python
import numpy as np
import pertpy as pt
import scanpy as sc
import pandas as pd


# Read data
adata = sc.read_h5ad('adata.h5ad')

# Initialize Milo analysis object
milo = pt.tl.Milo()
mdata = milo.load(adata)

# Build neighbor graph
sc.pp.neighbors(
    mdata["rna"],
    use_rep="X_umap", # Use existing UMAP dimensionality reduction
    n_pcs=20,         # Number of principal components, corresponds to d in miloR
    n_neighbors=30,   # Number of neighbors, corresponds to k in miloR
)

# Create neighborhoods
milo.make_nhoods(mdata["rna"], prop=0.1)
mdata = milo.count_nhoods(mdata, sample_col="sample_id")
# Group specification is more intuitive than in miloR, but order differs - Control comes first
mdata["rna"].obs["group"] = mdata["rna"].obs["group"].cat.reorder_categories(["Control", "Case"])

# Differential abundance analysis (using pydeseq2 solver)
milo.da_nhoods(mdata, design="~group", solver="pydeseq2")

# Build neighborhood graph for visualization
milo.build_nhood_graph(mdata)

# Filter significant results
significant_results = mdata['milo'].var.query('SpatialFDR < 0.1')

# All computations completed
```

## Postscript

In terms of actual performance, there isn't a huge difference between the two implementations. However, there are differences in computation results since the Python version is an algorithmic implementation rather than a step-by-step reproduction. Regarding computational consistency, I'll need to find a dataset from cellxgene to perform actual calculations in the future.
