---
title: 为了复现结果，需要一并复现bug
categories: Bioinformatics
date: 2025-11-26 21:11:53
tags: ['maftools', 'R']
---


只要活的够久，就会见到足够多的囧事情。
<!-- more -->

跟客户合作的一篇文章需要小修，我们需要根据审稿人的意见补充一些展示WES状况的图，由于之前负责的同事已经离职了，这项工作自然交给我了...

但是就像大部分半路接收的工作一样，拿到数据，一跑，结果不出所料的，跟文章里的图不一样。

需要绘制的是`maftools`产生的`oncoplot`，大概选择了20+个基因来展示配对样本的突变一致性，然而同样的文件，我画出来的突变类型跟前同事提供的基本都对不上。相比前同事的版本，我的图里`Splice_Site`消失了，同时虽然都检出突变，但是我的版本显示出大量的`MultiHit`，但是前同事的版本却基本都是每个基因单突变类型。

遂查看了原始文件，发现原始文件中确实有`Splice_Site`，所以判断是我的版本有问题，再次开始翻Github的Issue，经过一段折腾，没想到，涉及了2个bug，一个是我版本里的，另一个，是前同事的...

## 背景补充

首先补充一下，软件版本的差异。

前同事是使用个人电脑分析的，自行装的R和Rstudio，版本反正不可考了。我这边入职的是有已经购置了服务器，所以一直是使用服务器分析。然后因为我按照Pixi的原设计逻辑，一项目一配置，且我跨项目不保留`pixi.lock`，`pixi.toml`中除非有兼容性问题，并不限定版本，所以我用的包都会是依赖满足范围内，`bioconda`上有的最新版本。目前`maftools`在上面更新到`2.22`。

再然后是输入文件的情况，现在的单位使用的是ANNOVAR（跳槽了一圈又回到ANNOVAR了）

## `Splice_Site`解析问题
追溯我的版本中，`Splice_Site`的丢失情况，可以发现，是将ANNOVAR结果中的突变类型，映射成MAF对应的类型时出的问题，具体见下：

```R
ann[, `:=`(Hugo_Symbol, unlist(data.table::tstrsplit(Gene.refGene, 
    split = ";", keep = 1)))]
annovar_values = c(exonic = "RNA", splicing = "Splice_Site", 
    UTR5 = "5'UTR", UTR3 = "3'UTR", intronic = "Intron", 
    upstream = "5'Flank", downstream = "3'Flank", intergenic = "IGR", 
    `frameshift insertion` = "Frame_Shift_Ins", `frameshift deletion` = "Frame_Shift_Del", 
    `frameshift block substitution` = "Frameshift_INDEL", 
    `frameshift substitution` = "Frameshift_INDEL", stopgain = "Nonsense_Mutation", 
    stoploss = "Nonstop_Mutation", startloss = "Translation_Start_Site", 
    startgain = "Unknown", `nonframeshift insertion` = "In_Frame_Ins", 
    `nonframeshift deletion` = "In_Frame_Del", `nonframeshift block substitution` = "Inframe_INDEL", 
    `nonframeshift substitution` = "Inframe_INDEL", `nonsynonymous SNV` = "`Missense`_Mutation", 
    `synonymous SNV` = "Silent", unknown = "Unknown", ncRNA_exonic = "RNA", 
    ncRNA_intronic = "RNA", ncRNA_UTR3 = "RNA", ncRNA_UTR5 = "RNA", 
    ncRNA = "RNA", ncRNA_splicing = "RNA")
ann[, `:=`(Func.refGene, unlist(data.table::tstrsplit(x = as.character(Func.refGene), 
    split = ";", keep = 1)))]
ann[, `:=`(ExonicFunc.refGene, unlist(data.table::tstrsplit(x = as.character(ExonicFunc.refGene), 
    split = ";", keep = 1)))]
ann$Variant_Classification = ifelse(
    # 出问题的地方在这里，这里的逻辑是：
    # 当ExonicFunc.refGene不存在数值时，使用Func.refGene数值带入命名向量以进行映射，
    # 否则使用ExonicFunc.refGene带入映射。但是，根据ANNOVAR的文档，进行注释时一般都会
    # 加上 `-nasting .`，即用点号表示缺失，因此ExonicFunc.refGene根本不会缺失
    is.na(ann$ExonicFunc.refGene),  
    annovar_values[ann$Func.refGene],
    annovar_values[ann$ExonicFunc.refGene]
    )
```

经过上面描述的问题，由于ExonicFunc.refGene根本不会缺失，所有`ExonicFunc.refGene`为缺失的变异，其实也就是外显子区域以外的变异，都会因为`Variant_Classification`映射出缺失值被过滤掉，体现出来的效果就是，这些突变解析后全部不见了。

所以修复也很简单，额外判断一下`.`，或者前面`fread`读取的时候指定`.`是缺失就行了。

```R
ann$Variant_Classification = ifelse(
    # 我选择直接改判断
    is.na(ann$ExonicFunc.refGene) | ann$ExonicFunc.refGene == ".", 
    annovar_values[ann$Func.refGene],
    annovar_values[ann$ExonicFunc.refGene]
    )
```

这个问题追溯起来，其实不一定是maftools的原因，我往前倒了好多个版本，没见这个地方有修改，说不定是`fread`改了默认的`nastring`参数...反正我的问题解决了，我也不继续追了...

## `MultiHit`解析问题

这个问题是源自一个有[Issue](https://github.com/PoisonAlien/maftools/issues/347)的问题了。简单来说就是，maftools的老版本中，`MultiHit`的逻辑有问题。`MultiHit`想展示的，是一个基因发生两个不同的改变氨基酸的突变（`Missense` + `Missense` = MultiHit）。但是老版本中，实际上同一个基因发生不同类型的突变才会是`MultiHit`，如果都是`Missense`突变，则会直接显示`Missense`。

在我使用的最新版本中，这个问题已经被修复，因此我的结果中出现了大量的`MultiHit`...

那...能怎么办呢，对着Issue的修改，把Bug复现出来呗...

## 后记

前两天才吐槽过，`Azimuth`的Bug，我真没想到这么快，用过好多次的`maftools`也会有问题。

只能说敝人不只是天煞孤星，还是特么是扫把星...