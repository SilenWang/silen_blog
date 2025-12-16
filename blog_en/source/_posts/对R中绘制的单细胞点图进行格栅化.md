---
title: Rasterizing Single‑Cell Dot Plots in R
categories: Bioinformatics
date: 2025-12-17 01:08:35
tags: ['ggplot2', 'R', 'rasterization']
---

In bioinformatics visualization, we often need to handle plots containing tens of thousands of data points, such as scatter plots from single‑cell RNA‑seq data. When saved in vector formats like PDF, such graphics can suffer from huge file sizes and slow rendering (most software other than AI will simply freeze), because a vector file records the coordinates, color, size, and other attributes of every single point, resulting in a PDF with an enormous number of objects that hampers viewing and editing efficiency.

<!-- more -->

In contrast, bitmap images (e.g., PNG, JPEG) store the picture as a pixel matrix, which yields relatively fixed file sizes and fast rendering, but they lose quality when zoomed in. To combine the strengths of both approaches, we can **rasterize** the data‑dense layers of a plot—i.e., convert only specific layers to a bitmap while keeping other layers (such as axes, text labels) as vectors. This significantly reduces file size while preserving clarity where it matters most.

The following R code demonstrates how to use the `rasterise` function from the `ggrastr` package to rasterize layers of a `ggplot2`‑based graphic:

```R
library(ggrastr)
library(scRNAtoolVis)

# Create the plot; the same principle applies to any function that returns a ggplot2 object       
plot <- jjVolcano(
    diffData = markerspbmc.markers,
    log2FC.cutoff = 0.5, 
    col.type = "adjustP", 
    size  = 3.5,
    fontface = 'italic'
)

layers <- plot$layers
# After inspecting the layers, the third layer is the dot layer I want to rasterize
layers[[3]] <- rasterise(layers[[3]], dpi = 300)
# Replace the layers back into the plot
plot$layers <- layers
```

This technique is especially useful for high‑data‑volume visualizations such as single‑cell dot plots, genome‑browser tracks, etc. It retains the crispness of important annotations in vector form while making the saved file editable on ordinary computers (though if you intend to edit the points themselves, that’s another story…).
