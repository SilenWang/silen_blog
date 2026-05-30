---
title: Pitfalls When Plotting with Pathview
categories: Bioinformatics
date: 2026-05-28 16:00:00
tags: ['R', 'pathview', 'KEGG', 'data visualization', 'DEG analysis']
---

Recently, while doing differential gene expression analysis, I needed to map the expression values of DEGs onto KEGG pathways. I used Bioconductor's `pathview` package, but the resulting nodes were ridiculously large and the colors were wrong. With the help of AI, I fixed the issues and also picked up a new patching technique along the way.

<!-- more -->

## A Cleaner Way to Patch Functions in R

{% post_link 对发现bug的R函数进行热修复 [Previously, I fixed a buggy R function] %} by replacing the entire function, but when only a few lines need changing, that approach is neither elegant nor does it make the change obvious.

This time Deepseek used a different approach: `trace`.

`trace()` is an R debugging tool originally designed for tracing function calls, but it can also inject extra code into a function.

For example, to print `'123'` before `pathview`'s `render.kegg.node` executes:

```R
trace(pathview:::render.kegg.node,
      tracer = quote({
        print('123')
      }),
      print = FALSE)
```

If a function bug can be fixed by modifying certain settings before it runs, using `trace` is sufficient — no need to fully replace and re-register the function.

All three problems I encountered were solved this way.

## Case 1: Compound Node Size Explosion

I needed to plot the fatty acid metabolism pathway (ID 01212, a metabolic pathway, not a signaling pathway). The nodes are metabolites, and enzymes/genes are on the edges. `pathview` by default redraws metabolite nodes as circles and labels them with gene names, but the metabolite nodes were rendered way too large, obscuring edges and text.

Investigation revealed that `render.kegg.node` computes compound node dimensions incorrectly, producing excessive width and height values. Using `trace()` to intercept and fix:

```R
trace(pathview:::render.kegg.node,
      tracer = quote({
        if (type == "compound" && !is.null(plot.data)) {
          plot.data$width <- rep(2, nrow(plot.data))
          plot.data$height <- rep(2, nrow(plot.data))
        }
      }),
      print = FALSE)
```

How it works: `trace()` injects code before the function runs, forcing all metabolite node widths and heights to `2`.

## Case 2: Color Mapping Failure for Gene Boxes

Pathway 01212 is not a signaling pathway and has no gene nodes, so `pathview` cannot display up/down regulation using colored gene node boxes.

My idea was to encode the expression change in the gene label color instead. The `cols.ts` parameter in `render.kegg.node` is normally used to color certain elements. Deepseek set `text.col` to `cols.ts`, so text color reflects the logFC:

```R
if (type == "gene" && !is.null(cols.ts)) {
  text.col <- cols.ts
}
```

The resulting figure isn't pretty, but at least it's readable...

## Case 3: Different Pathways Need Different Rendering Strategies

Fixing metabolic pathways created a problem for signaling pathways. Signaling pathways have gene nodes that are normally rendered with color boxes. If the text color mapping also applies to signaling pathways, the text blends in with the background box and becomes invisible...

There was no reliable source for determining the pathway type, so I hardcoded it:

```R
text_color_pathways <- c("mcf01212", "mmu01212")
.pv_use_text_color <<- full_pathway_id %in% text_color_pathways
if (exists(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE) &&
    isTRUE(get(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE))) {
  text.col <- cols.ts
}
```

A global variable is used here to pass information to the callback injected by `trace` (because `trace`'s `quote` executes inside `pathview`'s namespace, so parameters cannot be passed directly).

## Other Parameter Tuning Notes

Several other parameters were also adjusted during plotting:

- `kegg.native = TRUE`: If set to `FALSE`, the pathway is rendered as a Graphviz directed graph that only preserves topology, which is basically unreadable when there are many nodes.
- `same.layer = FALSE`: When plotting a node-heavy pathway like 01212, `pathview` crashes unless this is set to `FALSE`.
- `limit = list(gene = max(abs(log2fc), na.rm=TRUE))`: Makes the color mapping range match the actual data; the default scales to -1 ~ 1.

## Acknowledgment

Thanks to the [discussion](https://www.biostars.org/p/385963/) on Biostars, which mentioned using `trace()` to modify `pathview` functions. Deepseek built on that to produce the solution.
