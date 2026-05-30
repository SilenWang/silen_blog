---
title: 使用Pathview绘图踩坑记录
categories: Bioinformatics
date: 2026-05-28 16:00:00
tags: ['R语言', 'pathview', 'KEGG', '数据可视化', 'DEG分析']
---

最近做差异基因分析，需要将差异基因表达量映射到 KEGG 通路上生成 pathway 图。我用 Bioconductor 的 `pathview` 包出图后发现节点大得离谱、颜色也不对。在 AI 帮助下解决了问题，还学到一种新的打补丁方式。

<!-- more -->

## R 中更简洁的打补丁方式

{% post_link 对发现bug的R函数进行热修复 [之前我修复 R 函数 Bug] %} 用的是整个替换函数的方式，但只需修改少数几行时，这种方式不够优雅，也不易看出改了什么。

这次 Deepseek 用了另一种方式：`trace`。

`trace()` 是 R 的调试工具，原本用于函数调用追踪，也可以往函数中插入额外代码。

例如，在 `pathview` 的 `render.kegg.node` 执行前打印 `'123'`：

```R
trace(pathview:::render.kegg.node,
      tracer = quote({
        print('123')
      }),
      print = FALSE)
```

如果函数 bug 可通过在运行前修改某些设置来修复，用 trace 即可，无需完整替换函数再注册。

我碰到的三个问题都靠这个方法解决。

## 用例一：化合物节点尺寸爆炸

我需要绘制编号 01212 的脂肪酸代谢通路（代谢通路，非信号通路），节点是代谢物，酶/基因在边上。`pathview` 默认将代谢物节点重绘为圆形并标记基因名称，但代谢物节点被渲染得巨大，遮住了边和文字。

排查发现 `render.kegg.node` 计算化合物节点尺寸时出错，导致宽高值过大。用 `trace()` 拦截修改：

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

原理：用 `trace()` 在函数前注入设置，将所有代谢物节点宽高强制设为 `2`。

## 用例二：基因盒子的颜色映射失效

01212 不是信号通路，没有基因节点，`pathview` 无法在基因节点上用颜色方框表示上下调。

我的思路是把颜色体现在基因文字颜色上。`render.kegg.node` 的 `cols.ts` 参数本用于给某些元素上色，Deepseek 将 `text.col` 设为 `cols.ts`，文字颜色就能反映 logFC 变化：

```R
if (type == "gene" && !is.null(cols.ts)) {
  text.col <- cols.ts
}
```

这样出来的图虽然不太好看，但好歹能看了...

## 用例三：不同通路需要不同的渲染策略

前面解决了代谢通路的问题，却给信号通路挖了坑。信号通路有基因节点，正常会渲染颜色，若文字颜色映射也对信号通路起效，文字和背景框颜色一致就看不到了...

判断通路类型没有好的信息来源，于是硬编码：

```R
text_color_pathways <- c("mcf01212", "mmu01212")
.pv_use_text_color <<- full_pathway_id %in% text_color_pathways
if (exists(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE) &&
    isTRUE(get(".pv_use_text_color", envir = .GlobalEnv, inherits = FALSE))) {
  text.col <- cols.ts
}
```

这里用全局变量传递信息给 `trace` 注入的回调（`trace` 的 `quote` 在 `pathview` 命名空间内执行，无法直接传参）。

## 其他调参记录

绘制中还调整了几个参数：

- `kegg.native = TRUE`：若用 `FALSE`，渲染成 Graphviz 有向图，只保留拓扑关系，节点多时基本不可读。
- `same.layer = FALSE`：绘制 01212 这种节点过多的图时，`pathview` 会崩溃，必须设为 `FALSE`。
- `limit = list(gene = max(abs(log2fc), na.rm=TRUE))`：让颜色映射范围与实际数据一致，默认缩放到 -1 ~ 1。

## 感谢

感谢 Biostars 上的[讨论](https://www.biostars.org/p/385963/)，提到了用 `trace()` 修改 `pathview` 函数的手法，Deepseek 据此写出了解决方案。