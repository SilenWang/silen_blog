---
title: 在R和Python下运行milo的对比
categories: Bioinformatics
date: 2025-10-28 15:57:51
tags: [单细胞分析, 差异丰度分析, Milo, 生物信息学]
---

Milo是一种用于单细胞RNA测序数据的差异丰度分析方法，能够检测不同条件下细胞邻域的组成变化。
<!-- more -->

## milo算法简介

[milo](https://www.nature.com/articles/s41587-021-01033-z)是一种专门为单细胞RNA测序数据设计的差异丰度分析方法。其核心思想是通过构建细胞邻域（neighborhoods）来检测不同实验条件下细胞群体组成的变化。

个人理解这个算法的设计目的是在不预先定义的情况下，找出两种情况下有群体比例差异的细胞，有了检测结果后，可以将这些有差异的细胞单独提出来进行进一步的分析。

## R环境下的miloR实现

R语言中使用miloR包进行分析，以下是相应的实现：

```r
# 加载所需库
library(Seurat)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(dplyr)
library(tidyr)

# 读取Seurat对象并转换为SingleCellExperiment
seurat_obj <- readRDS("seurat.rds")
sce <- as.SingleCellExperiment(seurat_obj)

# 设置随机种子保证结果可重复，本方法中主要是kNN部分有随机成分
set.seed(4466)

# 创建Milo对象
milo_obj <- Milo(sce)

# 构建领接结构，注意合理的k是kNN中的k
# 而d是dimension，即需要使用多少个降维维度/成分
milo_obj <- buildGraph(
    milo_obj,
    k = 30, d = 20,
    transposed = TRUE, 
    reduced.dim = "UMAP"
)

# 定义邻域，注意这里的k和d以及reduced_dims都要跟前一步一致
# 这里的refinement_scheme是为了配合后续使用一种快速算法
# 在数据集中细胞非常多时，有明显的计算加速作用
milo_obj <- makeNhoods(
    x = milo_obj,
    prop = 0.1, k = 30, d = 20,
    refined = TRUE,
    reduced_dims = "UMAP",
    refinement_scheme = "graph"
)

# 计数邻域中的细胞
milo_obj <- countCells(
    milo_obj,
    meta.data = data.frame(colData(milo_obj)),
    sample = "sample_id"
)

# 准备实验设计，其实就是构建一个有样品编号和分组tag的数据框
traj_design <- data.frame(colData(milo_obj))[,c("sample_id", "group")] %>%
    distinct() %>%
    mutate(
        sample_id = as.character(sample_id),
        group = as.character(group)
    )
rownames(traj_design) <- traj_design$sample_id
traj_design <- traj_design[colnames(nhoodCounts(milo_obj)), , drop=FALSE]

# 差异丰度检验，这里的fdr.weighting = "graph-overlap"
# 配合前面makeNhoods中的参数，提高计算速度
# 注意这里可以进行比较组的指定，指定的方式比较特殊，是用类公式的语法进行指定
# design = ~ 0 + group 代表使用group变量进行分组，前面的 0 + 是固定写法
# model.contrasts = c('groupCase - groupControl') 的意义是
# 使用 group 列等于 Case 的行作为目标，等于 Control的行作为对照
da_results <- testNhoods(
    milo_obj,
    design = ~ 0 + group,
    model.contrasts = c('groupCase - groupControl'),
    design.df = traj_design,
    fdr.weighting = "graph-overlap"
)

# 查看显著结果
da_results %>%
  arrange(SpatialFDR) %>%
  filter(SpatialFDR < 0.1) %>%
  tail(3)

# 构建邻域图
milo_obj <- buildNhoodGraph(milo_obj)

# 所有计算步骤完成
# 后续可使用plotUMAP或其他函数进行可视化展示
```

## Python环境下的Milo实现

Python中有多个milo实现，除了官方[miloR项目](https://github.com/MarioniLab/miloR)中提到的[Mipopy](https://github.com/emdann/milopy)，[pertpy](https://github.com/scverse/pertpy)库也实现Milo的算法，以下是用pertpy完成计算的例子

```python
import numpy as np
import pertpy as pt
import scanpy as sc
import pandas as pd


# 读取数据
adata = sc.read_h5ad('adata.h5ad')

# 初始化Milo分析对象
milo = pt.tl.Milo()
mdata = milo.load(adata)

# 构建邻居图
sc.pp.neighbors(
    mdata["rna"],
    use_rep="X_umap", # 利用已有umap降维信息
    n_pcs=20,         # 主成分数，对应miloR的d参数
    n_neighbors=30,   # 邻居数，对应miloR的k参数
)

# 创建邻域
milo.make_nhoods(mdata["rna"], prop=0.1)
mdata = milo.count_nhoods(mdata, sample_col="sample_id")
# 这里指定Case Control的方式比miloR更容易懂一些，不过顺序有些差异，这里Control在前
mdata["rna"].obs["group"] = mdata["rna"].obs["group"].cat.reorder_categories(["Control", "Case"])

# 差异丰度分析（使用pydeseq2求解器）
milo.da_nhoods(mdata, design="~group", solver="pydeseq2")

# 构建邻域图，用于可视化
milo.build_nhood_graph(mdata)

# 筛选显著结果
significant_results = mdata['milo'].var.query('SpatialFDR < 0.1')

# 完成所有计算
```

## 后记

实际运行的性能上，两者差距木有太大，但是计算结果上还是又差的，毕竟py版本是算法实现，而不是一点点复现。具体的计算一致性几何，后续我从cellxgene上找个数据集再来实际算算好了。