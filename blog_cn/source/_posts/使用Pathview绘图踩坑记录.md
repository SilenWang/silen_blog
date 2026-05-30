---
title: 使用Pathview绘图踩坑记录
categories: Bioinformatics
date: 2026-05-28 16:00:00
tags: ['R语言', 'pathview', 'KEGG', '数据可视化', 'DEG分析']
---

最近在跑差异基因分析的流程，有个环节是把差异基因的表达量映射到 KEGG 通路上，生成那种直观的 pathway 图。我用的是 Bioconductor 上的 `pathview` 包，结果出图一看，好家伙，有些节点大得离谱，而且颜色也不对劲。

<!-- more -->

## 问题一：化合物节点尺寸爆炸

第一次跑 `pathview`，出来的图长这样——某些化合物（Compound）节点被渲染得巨大无比，严重挤压了其他节点，整张图根本没法看。

查了一下，问题是 `pathview` 内部的 `render.kegg.node` 函数在计算化合物节点尺寸时出了问题，导致宽度和高度被设成了很大的值。这个函数不是导出的 API，但可以通过 `trace()` 来拦截和修改它的行为：

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

原理很简单：用 `trace()` 注入一个回调，在 `render.kegg.node` 执行到我们插入的代码位置时，将所有化合物节点的宽高强制设为 `2`。

## 问题二：基因盒子的颜色映射失效

解决尺寸问题后，新的问题又来了。对于某些通路，比如 `mmu01212`（脂肪酸代谢通路），`pathview` 默认会在基因节点上画带颜色的方框来表示上下调——但其实这个通路在 KEGG 的 XML 里并没有标准的 "gene box" 结构，所以颜色根本不会显示。

那怎么办呢？我的思路是把颜色的信息体现在基因名字的文字颜色上。`render.kegg.node` 里有一个 `cols.ts` 参数，它本身就是为了给某些元素上色用的。只要把 `text.col` 设成 `cols.ts`，文字颜色就能正确反映 logFC 的变化了：

```R
if (type == "gene" && !is.null(cols.ts)) {
  text.col <- cols.ts
}
```

注意这只是配色方案的一种——有的通路有标准的 gene box，用默认的方框上色效果更好；有的通路（如 `mmu01212`、`mcf01212`）没有，就得用文字颜色。

## 问题三：不同通路需要不同的渲染策略

这就引出了第三个问题：不同通路需要不同的渲染方式。我总不能一刀切——让所有的通路都用文字颜色，那样反而会让好端端的 gene box 配色的通路变得奇怪。

解决方案是定义一个白名单列表，只对特定的通路采用文字颜色方案：

```R
text_color_pathways <- c("mcf01212", "mmu01212")
```

然后在通路渲染循环里，根据当前通路是否在白名单中来决定是否启用文字颜色：

```R
.pv_use_text_color <<- full_pathway_id %in% text_color_pathways
```

注意这里用了一个全局变量来把信息传递给 `trace` 注入的回调函数（因为 `trace` 里的 `quote` 是在 `pathview` 的命名空间内执行的，没法直接传参）：

```R
if (exists(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE) &&
    isTRUE(get(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE))) {
  text.col <- cols.ts
}
```

## 彩蛋：标记显著差异基因

既然都已经在上一步拿到了所有基因的差异表达数据，那不如顺便把显著差异基因（padj < 0.05）也标出来。我在渲染循环里又加了一句，把 padj 信息存到全局变量里：

```R
.pv_gene_padj <<- setNames(deg_data$padj, as.character(deg_data$gene_id))
```

然后在 `trace` 回调里，对于显著差异的基因，在基因名后面加个星号 `*`：

```R
if (exists(".pv_gene_padj", envir = .GlobalEnv, inherits = FALSE)) {
  padj_vec <- get(".pv_gene_padj", envir = .GlobalEnv, inherits = FALSE)
  if (!is.null(padj_vec) && length(padj_vec) > 0) {
    node_ids <- gsub("^[a-z]+:", "", as.character(plot.data$kegg.names))
    sig <- node_ids %in% names(padj_vec) & padj_vec[node_ids] < 0.05
    sig[is.na(sig)] <- FALSE
    plot.data$labels[sig] <- paste0(plot.data$labels[sig], "*")
  }
}
```

这样，一眼就能看出哪些基因是显著差异的，非常直观。

## 其他调参记录

除了上面三个大坑，还调整了几个参数：

- `kegg.native = TRUE`：之前用的 `FALSE`，会渲染成 Graphviz 布局的图，看着不太对劲。改成 `TRUE` 使用 KEGG 原生的 PNG 底图，然后把数据覆盖上去，效果好了很多。
- `map.symbol = FALSE`：如果是 `TRUE` 会把 Entrez ID 转成 gene symbol，但很多时候并不需要这个转换。
- `same.layer = TRUE`：让基因数据和通路底图在同一层渲染，避免分层导致的错位问题。
- `limit = list(gene = max(abs(log2fc), na.rm=TRUE))`：让颜色映射的范围跟数据实际范围一致，而不是用默认值，这样颜色对比更明显。

## 后记

这些修改都记录在 `VyBioTx/BioWorkflow` 仓库的三个 commit 里：

- `ee26122`：修复化合物节点尺寸
- `7542a82`：增加文字颜色表示 logFC
- `68ffddf`：支持逐通路配置 + 添加显著性标记

还要感谢 Biostars 上的这个讨论（<https://www.biostars.org/p/385963/>），里面提到了用 `trace()` 来修改 `pathview` 函数的手法，给了我很大启发。

说实话，这种用 `trace()` 做热修复的方式其实挺 hacky 的，但毕竟不能指望 R 包作者满足所有人的需求。在不修改源代码的前提下，这种方式可以灵活地扩展包的渲染能力，对数据分析流程来说已经够用了。
