---
title: 使用GNU Parallel在Bash中执行并行
categories: Script
date: 2019-03-31 23:07:26
tags: ["bash", "并行计算"]
---

之前已经在R和Python中都使用过并行了, 我是最近才知道原来Bash下面也有简单好用的并行程序: *GNU Parallel*[^1]. 这是一个Perl实现的程序, 可以方便的将Bash中各种命令并行执行.

[^1]: https://doi.org/10.5281/zenodo.1146014

<!-- 摘要部分 -->
<!-- more -->

由于刚刚上手, 现在暂时只知道一种使用方式, 就是替代For循环, 将并行的命令替换成并行执行, 比如下面写了一个分染色体合并vcf的例子:

```bash
for chr in {{1..22},X,Y};do
    bcftools view -f PASS -r chr${chr} sample.vcf > sample.chr${chr}.flt.vcf
done
bcftools merge sample.chr{{1..22},X,Y}.vcf > sample.flt.vcf
```

使用GNU Parallel替代可以是这样子:

```bash
echo {{1..22},X,Y} | tr " " "\n" | parallel 'bcftools view -f PASS -r chr{} sample.vcf > sample.chr{}.flt.vcf'
bcftools merge sample.chr{{1..22},X,Y}.vcf > sample.flt.vcf
```

上面的会是循环`bcftools view`然后合并, 下面则是并行执行后再合并. 当然我这个例子没有意义...因为第一种先循环再合并...和直接过滤没有什么区别...

默认情况下`parallel`有多少核心就用多少核心, 所以实际使用可能会使用`-j`指定核心上限或者使用`--load 80%`指定用80%的核心.

有了这个东西, 只要使用的工具是支持输入参数定位分析区域的, 都可以快速的将其变成并行执行, 加快分析速度. 这么一来, 用Snakemake写流程舒心了很多...

关于这个程序的其他用法之后有缘继续更...