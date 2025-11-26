---
title: To Reproduce Results, I Had to Reproduce Bugs Too
categories: Bioinformatics
date: 2025-11-26 21:11:53
tags: ['maftools', 'R']
---

If you live long enough, you'll encounter plenty of awkward situations.
<!-- more -->

An article we're collaborating on with a client needed minor revisions, and I had to add some plots showing WES status according to the reviewers' comments. Since the colleague previously responsible had left the company, this task naturally fell to me...

But like most work taken over midway, when I got the data and ran it, unsurprisingly, the results didn't match the figures in the article.

We needed to generate `oncoplot` using `maftools`, selecting about 20+ genes to show mutation consistency in paired samples. However, with the same input files, the mutation types in my plots were fundamentally different from those provided by my former colleague. Compared to her version, `Splice_Site` mutations disappeared in my plots, and while mutations were detected, my version showed lots of `MultiHit`, whereas her version mostly showed single mutation types per gene.

I checked the original files and found that `Splice_Site` mutations were indeed present, so I concluded the issue was with my version. I started digging through GitHub issues, and after some investigation, I discovered it involved two bugs - one in my version, and another in my former colleague's...

## Background Information

First, let me explain the software version differences.

My former colleague used her personal computer for analysis, with self-installed R and RStudio - the versions are unknown. I work on a company-purchased server, so I've been using server-based analysis. Following Pixi's original design logic, I use per-project configurations and don't preserve `pixi.lock` across projects. In `pixi.toml`, unless there are compatibility issues, I don't specify versions, so the packages I use are the latest available on `bioconda` within dependency constraints. Currently, `maftools` is updated to `2.22` there.

Regarding input files, our current workplace uses ANNOVAR (I've come full circle back to ANNOVAR after job hopping).

## The `Splice_Site` Parsing Issue
Tracing the disappearance of `Splice_Site` in my version revealed the problem occurred when mapping mutation types from ANNOVAR results to MAF-compatible types, as shown below:

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
    # The problem is here. The logic is:
    # When ExonicFunc.refGene has no value, use Func.refGene value for mapping via named vector,
    # Otherwise use ExonicFunc.refGene for mapping. However, according to ANNOVAR documentation,
    # annotations typically use `-nasting .`, meaning dots represent missing values,
    # so ExonicFunc.refGene would never actually be missing
    is.na(ann$ExonicFunc.refGene),  
    annovar_values[ann$Func.refGene],
    annovar_values[ann$ExonicFunc.refGene]
    )
```

Due to the described issue, since ExonicFunc.refGene would never actually be missing, all mutations where `ExonicFunc.refGene` appeared missing (actually mutations outside exon regions) would be filtered out because `Variant_Classification` would map to missing values. The effect was that these mutations completely disappeared after parsing.

The fix was simple - either check for `.` explicitly, or specify during `fread` that `.` represents missing values.

```R
ann$Variant_Classification = ifelse(
    # I chose to modify the condition directly
    is.na(ann$ExonicFunc.refGene) | ann$ExonicFunc.refGene == ".", 
    annovar_values[ann$Func.refGene],
    annovar_values[ann$ExonicFunc.refGene]
    )
```

Tracing this issue, it might not necessarily be maftools' fault. I checked several previous versions and didn't see changes in this part - maybe `fread` changed the default `nastring` parameter... Anyway, my problem was solved, and I didn't investigate further.

## The `MultiHit` Parsing Issue

This issue originates from a [reported GitHub issue](https://github.com/PoisonAlien/maftools/issues/347). Simply put, the logic for `MultiHit` in older versions of maftools was flawed. `MultiHit` was intended to show when a gene had two different amino acid-changing mutations (`Missense` + `Missense` = `MultiHit`). However, in older versions, only different types of mutations in the same gene would be classified as `MultiHit`. If both were `Missense` mutations, it would simply display as `Missense`.

In the latest version I was using, this issue had been fixed, which is why my results showed lots of `MultiHit`...

So... what could I do? I had to reproduce the bug according to the issue's description...

## Afterword

I had just complained about `Azimuth`'s bugs a couple of days ago, and I never expected that `maftools`, which I've used many times, would also have issues.

I guess I'm not just an unlucky star, but also a freaking jinx...