---
title: 绘制一个配对样本的Circos图
categories: Bioinformatics
date: 2025-10-06 22:09:54
tags: [R, circlize, 基因组, 突变, CNV]
---

在癌症研究中，比较肿瘤样本和类器官模型的基因组特征对于验证模型的可靠性至关重要。Circos图能比较直观地展示基因组中检测到突变的情况，所以在描述代表样本的总体检测结果时非常常用。在阅读circlize的文档的时候，看到作者给了一个配对样品的展示例子，我觉得用来展示配对的原代样本和类器官挺合适的，就搓了一个用来展示配对样本的图。代码主要参考自官方文档的[9.5 Concatenating two genomes](https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html#concatenating-two-genomes)

<!-- more -->

## 准备工作

### 加载必要的R包

首先，我们需要加载一些必要的R包：

```r
library(circlize)
library(dplyr)
library(stringr)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)
```

其中，`MutationalPatterns`主要是用来加载vcf文件，这样不用手动处理，`BSgenome.Hsapiens.UCSC.hg38`是绘制基本的染色体区间时需要包提供的信息

### 数据准备

我这次展示了三方面的信息，分别是：
- 体细胞突变检测结果（VCF文件）
- CNV检测结果（来自CNVKit）
- 测序覆盖度数据（直接利用CNV软件给的结果，这样不用重复）

## 构建基因组框架数据

我们首先构建一个包含肿瘤和类器官基因组的合并框架：

```r
# 去除XY，构建合并的基因组数据
tumor_cytoband <- read.cytoband(species = "hg38")$df %>% filter(!(V1 %in% c('chrX', 'chrY')))
organoid_cytoband <- read.cytoband(species = "hg38")$df %>% filter(!(V1 %in% c('chrX', 'chrY')))

tumor_cytoband[ ,1] <- paste0("tumor_", tumor_cytoband[, 1])
organoid_cytoband[ ,1] <- paste0("organoid_", organoid_cytoband[, 1])
cytoband <- rbind(tumor_cytoband, organoid_cytoband)
```

这里实际上是产生了名为cytoband的数据框，这个数据框丽包含要绘制的两个基因组的基本信息。

## 处理突变数据

使用`MutationalPatterns`包读取VCF文件，并转换为适合绘制的格式：

```r
vcfs <- read_vcfs_as_granges(
    c(
        'tumor.vcf',
        'pdo.vcf'
    ),
    c("tumor", "pdo"),
    "BSgenome.Hsapiens.UCSC.hg38"
)

mut_data <- rbind(
    data.frame(vcfs[[t_sample]]) %>%
        mutate(value = 1) %>%
        mutate(seqnames = str_c("tumor_", seqnames)) %>%
        select(seqnames, start, end, value),
    data.frame(vcfs[[o_sample]]) %>%
        mutate(value = 1) %>%
        mutate(seqnames = str_c("organoid_", seqnames)) %>%
        select(seqnames, start, end, value)
)
```

## 处理CNV和覆盖度数据

CNV和覆盖度的数据都取自CNVkit：

```r
# CNV数据
cnv_data <- rbind(
    read.table(
            "tumor.cnvkit.call.cns",
            sep="\t", header = F,
            col.names = c("chromosome", "start", "end", "cn")
        ) %>%
        mutate(value = cn) %>%
        mutate(chromosome = str_c("tumor_", chromosome)) %>%
        select(chromosome, start, end, value),
    read.table(
            "pdo.cnvkit.call.cns",
            sep="\t", header = F,
            col.names = c("chromosome", "start", "end", "cn")
        ) %>%
        mutate(value = cn) %>%
        mutate(chromosome = str_c("organoid_", chromosome)) %>%
        select(chromosome, start, end, value)
) %>% mutate(value = pmin(value, 4)) # 限制CNV值不超过4

# 覆盖度数据
cov_data <- rbind(
    read.table("tumor.cnvkit.cov.cnn", sep="\t", header = T) %>%
        mutate(value = log2) %>%
        mutate(chromosome = str_c("tumor_", chromosome)) %>%
        select(chromosome, start, end, value),
    read.table("pdo.cnvkit.cov.cnn", sep="\t", header = T) %>%
        mutate(value = log2) %>%
        mutate(chromosome = str_c("organoid_", chromosome)) %>%
        select(chromosome, start, end, value)
)
```

## 绘制Circos图

### 设置颜色和初始化

```r
# 颜色设置
red <- "#FFC6D6"
blue <- "#bebcff"
green <- "#9CCF83"
orange <- "#fbd988ff"
purpule <- "#ff726dff"
black <- "#000000"

# 设置图形大小
png(filename = "demo.circos.png", width = 1200, height = 1200, res = 200)

# 构建基础的circos框架
chromosome.index = c(
    paste0("tumor_chr", c(1:22)), 
    rev(paste0("organoid_chr", c(1:22)))
)

circos.par(gap.after = c(rep(1, 22), rep(1, 22)))
circos.initializeWithIdeogram(
    cytoband,
    plotType = NULL, 
    chromosome.index = chromosome.index
)
```

### 添加染色体轨道

```r
# 设置染色体track
circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
            gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
        },
    track.height = mm_h(1),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA
)
highlight.chromosome(paste0("tumor_chr", c(1:22)), col = red, track.index = 1)
highlight.chromosome(paste0("organoid_chr", c(1:22)), col = blue, track.index = 1)
circos.genomicIdeogram(cytoband)
```

### 添加覆盖度、CNV和突变密度轨道

```r
# 测序覆盖度
circos.genomicDensity(cov_data, col=orange, track.height = 0.1, window.size = 1e7)

# CNV轨道
circos.genomicTrackPlotRegion(
    cnv_data, 
    ylim = c(0, 4),
    panel.fun = function(region, value, ...) {
        cell.xlim = get.cell.meta.data("cell.xlim")
        for(h in c(0, 1, 2, 3, 4)) {
            circos.lines(cell.xlim, c(h, h), col = black)
        }
        col = ifelse(value[[1]] > 2, "red",
            ifelse(value[[1]] == 2, "green", "blue")
        )
        i = getI(...)
        circos.genomicRect(region, value, col = col, ytop = value + 0.3, ybottom = value - 0.3 , border = NA)
    },
    track.height = 0.1
)

# 突变密度
circos.genomicDensity(mut_data, col=purpule, track.height = 0.1, window.size = 1e7)
```

### 添加文本和图例

```r
# 设置中间的文字
text(0, 0.2, "Demo", cex = 2, font = 2)
# 设置两角的文字
text(-0.9, -0.8, "Tumor\nGenome")
text(0.9, 0.8, "Organoid\nGenome")

# 设置中间文字下方的图例
legend(
  x = 0,
  y = 0,
  legend = c("Coverage", "CNV gain", "CNV neutral", "CNV loss",  "Mutation density"),
  col = c(orange, "red", "green", "blue", purpule),
  pch = 15,
  pt.cex = 1,
  cex = 0.8,
  bty = "n",
  xjust = 0.5,
)

circos.clear()
dev.off()
```

## 完整代码

```r
library(circlize)
library(dplyr)
library(stringr)
library(MutationalPatterns)
library(BSgenome.Hsapiens.UCSC.hg38)

tumor_cytoband <- read.cytoband(species = "hg38")$df %>% filter(!(V1 %in% c('chrX', 'chrY')))
organoid_cytoband <- read.cytoband(species = "hg38")$df %>% filter(!(V1 %in% c('chrX', 'chrY')))

tumor_cytoband[ ,1] <- paste0("tumor_", tumor_cytoband[, 1])
organoid_cytoband[ ,1] <- paste0("organoid_", organoid_cytoband[, 1])
cytoband <- rbind(tumor_cytoband, organoid_cytoband)

t_sample <- "tumor"
o_sample <- "pdo"

vcfs <- read_vcfs_as_granges(
    c(
        'tumor.vcf',
        'pdo.vcf'
    ),
    c(t_sample, o_sample),
    "BSgenome.Hsapiens.UCSC.hg38"
)

mut_data <- rbind(
    data.frame(vcfs[[t_sample]]) %>%
        mutate(value = 1) %>%
        mutate(seqnames = str_c("tumor_", seqnames)) %>%
        select(seqnames, start, end, value),
    data.frame(vcfs[[o_sample]]) %>%
        mutate(value = 1) %>%
        mutate(seqnames = str_c("organoid_", seqnames)) %>%
        select(seqnames, start, end, value)
)


cnv_data <- rbind(
    read.table(
            "tumor.cnvkit.call.cns",
            sep="\t", header = T
        ) %>%
        mutate(value = cn) %>%
        mutate(chromosome = str_c("tumor_", chromosome)) %>%
        select(chromosome, start, end, value),
    read.table(
            "pdo.cnvkit.call.cns",
            sep="\t", header = T
        ) %>%
        mutate(value = cn) %>%
        mutate(chromosome = str_c("organoid_", chromosome)) %>%
        select(chromosome, start, end, value)
) %>% mutate(value = pmin(value, 4))


cov_data <- rbind(
    read.table("tumor.cnvkit.cov.cnn", sep="\t", header = T) %>%
        mutate(value = log2) %>%
        mutate(chromosome = str_c("tumor_", chromosome)) %>%
        select(chromosome, start, end, value),
    read.table("pdo.cnvkit.cov.cnn", sep="\t", header = T) %>%
        mutate(value = log2) %>%
        mutate(chromosome = str_c("organoid_", chromosome)) %>%
        select(chromosome, start, end, value)
)


########################################### 绘图 ###########################################

red <- "#FFC6D6"
blue <- "#bebcff"
green <- "#9CCF83"
orange <- "#fbd988ff"
purpule <- "#ff726dff"
black <- "#000000"

png(filename = "demo.circos.png", width = 1200, height = 1200, res = 200)

chromosome.index = c(
    paste0("tumor_chr", c(1:22)), 
    rev(paste0("organoid_chr", c(1:22)))
)

circos.par(gap.after = c(rep(1, 22), rep(1, 22))) # 设置区隔？
circos.initializeWithIdeogram(
    cytoband,
    plotType = NULL, 
    chromosome.index = chromosome.index
)

circos.track(
    ylim = c(0, 1),
    panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[2] + mm_y(2), 
            gsub(".*chr", "", CELL_META$sector.index), cex = 0.6, niceFacing = TRUE)
        },
    track.height = mm_h(1),
    cell.padding = c(0, 0, 0, 0),
    bg.border = NA
)
highlight.chromosome(paste0("tumor_chr", c(1:22)), col = red, track.index = 1)
highlight.chromosome(paste0("organoid_chr", c(1:22)), col = blue, track.index = 1)
circos.genomicIdeogram(cytoband)

circos.genomicDensity(cov_data, col=orange, track.height = 0.1, window.size = 1e7)

circos.genomicTrackPlotRegion(
    cnv_data, 
    ylim = c(0, 4),
    panel.fun = function(region, value, ...) {
        cell.xlim = get.cell.meta.data("cell.xlim")
        for(h in c(0, 1, 2, 3, 4)) { 
            circos.lines(cell.xlim, c(h, h), col = black)
        }
        col = ifelse(value[[1]] > 2, "red",
            ifelse(value[[1]] == 2, "green", "blue")
        ) 
        i = getI(...)
        circos.genomicRect(region, value, col = col, ytop = value + 0.3, ybottom = value - 0.3 , border = NA)
    },
    track.height = 0.1
)

circos.genomicDensity(mut_data, col=purpule, track.height = 0.1, window.size = 1e7)

text(0, 0.2, "Demo", cex = 2, font = 2)
text(-0.9, -0.8, "Tumor\nGenome")
text(0.9, 0.8, "Organoid\nGenome")

legend(
  x = 0,           
  y = 0,        
  legend = c("Coverage", "CNV gain", "CNV neutral", "CNV loss",  "Mutation density"),
  col = c(orange, "red", "green", "blue", purpule),
  pch = 15,        
  pt.cex = 1,
  cex = 0.8,
  bty = "n",       
  xjust = 0.5,     
)

circos.clear()

dev.off()
```

## 成品展示

- 这里使用的是一个样本生成的数据画了上下两圈，所以两边检测结果看上去是一样的，实际要是配对结果能这么一致，那高低得给佛祖多磕几个了。

![circos_demo](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/10/upgit_20251007_1759767842.jpg)