---
title: 修复pixi部署部分bioconda_r包后出现的缺依赖问题
categories: Bioinformatics
date: 2025-10-29 03:37:56
tags: [pixi, R, Bioconductor, 依赖管理]
---

在使用pixi管理生物信息学分析环境时，经常会遇到一些Bioconductor的R包安装后出现依赖缺失的问题。目前暂不清楚这个问题的原因，用了pixi一年了，这个问题到目前为止（2025.10）也木有修复，因此本文介绍如何通过pixi tasks功能来解决这类问题。

<!-- more -->

## 问题描述

在最近一年的pixi使用中，发现部署生物信息学分析环境时，`Seurat`和`maftools`这样的依赖`GenomeInfoDbData`、`BSgenome.Hsapiens.UCSC.hg38`获取基因组数据R包的依赖总是有问题，即使在依赖处显式的指定并安装这些包到pixi环境中，实际载入的时候还是会显示这些包并不存在，因此还需要配置额外的内容来补充这些依赖。

## 解决方案

我目前是通过pixi的tasks功能，将安装依赖的命令变成一个任务，以在部署环境后快速修复依赖问题，实际的配置文件示例如下：

```toml
[workspace]
name = "azimuth_demo"
version = "1.1"
description = "workspace with azimuth and jupyter"
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[environments]
label = ['kernel', 'label']

[feature.kernel.dependencies]
r-irkernel = '*'
jupyterlab = '*'

[feature.label.dependencies]
r-base = '*'
r-azimuth = '*'
r-BiocManager = '*'

[feature.label.tasks]
GenomeInfoDbData = {cmd = 'Rscript -e "BiocManager::install(\"GenomeInfoDbData\")"'}
BSgenome = {cmd = 'Rscript -e "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")"'}
EnsDb = {cmd = 'Rscript -e "BiocManager::install(\"EnsDb.Hsapiens.v86\")"'}
JASPAR2020 = {cmd = 'Rscript -e "BiocManager::install(\"JASPAR2020\")"'}
r_dep = {cmd = 'echo "bio dep for R done"', depends-on=['GenomeInfoDbData', 'BSgenome', 'EnsDb', 'JASPAR2020']}
```

激活pixi环境后，执行`pixi run r_dep`就可以完成依赖修复，然后开始工作了。
