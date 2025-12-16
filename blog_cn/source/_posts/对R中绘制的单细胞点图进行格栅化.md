---
title: 对R中绘制的单细胞点图进行格栅化
categories: Bioinformatics
date: 2025-12-17 01:08:35
tags: ['ggplot2', 'R', '格栅化', 'rasterization']
---

在生物信息学绘图中，我们经常需要处理包含成千上万个数据点的图形，例如单细胞RNA测序的散点图。这类图形在保存为PDF等矢量格式时会面临文件过大、渲染缓慢的问题（除了AI其他软件基本都会直接死机），因为矢量图会记录每个数据点的坐标、颜色、大小等属性，导致PDF文件包含大量对象，进而影响查看和编辑的效率。

<!-- more -->

相比之下，位图（如PNG、JPEG）将图像存储为像素矩阵，文件大小相对固定且渲染快速，但放大时会失真。结合两种格式的优点，我们可以对图形中数据密集的图层进行**格栅化（rasterization）**，即只将特定图层转换为位图，同时保持其他图层（如坐标轴、文字标签）为矢量，从而在保证清晰度的前提下显著降低文件大小。

下面这段R代码展示了如何使用`ggrastr`包的`rasterise`函数对`ggplot2`体系的图形对象进行格栅化：

```R
library(ggrastr)
library(scRNAtoolVis)

# 绘制图形，任何会返回ggplot2对象的绘图函数原理一样       
plot <- jjVolcano(
    diffData = markerspbmc.markers,
    log2FC.cutoff = 0.5, 
    col.type = "adjustP", 
    size  = 3.5,
    fontface = 'italic'
)

layers <- plot$layers
# 查看layers内容后，第三个图层是我想格栅化的点图
layers[[3]] <- rasterise(layers[[3]], dpi = 300)
# 替换回 plot
plot$layers <- layers
```

这种方法特别适合用于单细胞点图、基因组浏览器图等大数据量的可视化场景，既能保留重要标注的矢量清晰度，又能让存下来的图形可被一般电脑编辑（当然如果你要编辑点本身、当我没说...
