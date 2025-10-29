---
title: 解决scanpy读写loom文件格式丢失obs和var的index的问题
categories: Bioinformatics
date: 2025-10-28 17:22:12
tags: ['scanpy', 'loom', 'adata']
---

在我们现在基于scanp的单细胞流程中，有一步需要将AnnData对象保存为loom格式。但是与保存为h5ad不同，当我们不做任何处理，将AnnData对象写入loom文件后再次读取时，会发现`obs`和`var`的索引（index）信息丢失了，这些索引（通常是细胞条形码和基因名）变成了普通的数字编号。

<!-- more -->

## 问题描述

根据实际的测试和[write_loom方法的文档](https://anndata.readthedocs.io/en/latest/generated/anndata.AnnData.write_loom.html)，中的描述，可以推测实际的读写，`write_loom` 写入时，实际上的不会写入obs和var两个表中的index部分，即使通过参数指定`write_obsm_varm=True`，写入时也只会按照固定的名称分别产生两个新列）（obs_names和var_names）。

而根据[read_loom方法的文档](https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_loom.html)，`read_loom`在读取的时候，实际上要在obs中找名为`CellID`的列，在var中找名为`Gene`的列来分别作为两张表的编号。

这样一来一去，默认不存index，即使存了index，`write_loom`存的和`read_loom`读的也并不一致，就像两个部门对接没有做好，关键信息就这么没了...

## 解决方案

我的目的是尽可能不对原流程代码进行变更，因此直接产生原流程可读的loom文件是最好的方式：

```python
# 保存索引信息到obs和var的列中
adata.obs['CellID'] = adata.obs.index.tolist()
adata.var['Gene'] = adata.var.index.tolist()

# 写入loom文件（包含obsm和varm信息）
merge_adata.write_loom('/path/to/raw.loom')

# 读取loom文件，由于CellID和Gene存在，
# 读取的adata能成功保留两个表的index信息
adata = sc.read_loom('/path/to/raw.loom')
```
