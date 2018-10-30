---
title: snakemake的使用
categories: Bioinfomatic
date: 2018-10-30 12:40:03
tags: ['python', 'snakemake']
---

写流程是生物信息分析中经常需要涉及的工作, 一个成熟完善的好工具可以使工作的效率大大提高.
<!-- more -->

最简单的流程编写方式可以是将所有需要运行的命令都写成一个shell脚本, 并使用`&&`或其他方式指定任务的先后顺序, 然后挂到后台就好.

这么做虽然容易, 但是其重复性和抗错误性实在堪忧. 因此也就有了各种各样的流程管理工具: 我过去单位用的`sjm`, 10x Genomics自己开发的`Martian`, Broad Institute的`WDL`, 我最近在用的`snakemake`等等.

## snakemake简介
`snakemake`是德国xxx大学xxx实验室开发的生物信息分析流程管理工具, 其全部代码完全由python实现, 因此可以很好的和python程序整合在一起(当然其也有一套完整的命令行程序).

`snakemake`的特点包括但不限于:

- 简单易用
- 可直接使用python代码
- 内置分布式运算集群支持(SGE, LSF, etc.)
- 内置Conda环境部署支持
- 内置容器支持(注意这个容器不是docker)
- 自动处理依赖并进行断点重启
- 内置测试功能(benchmark)

## snakemake简单示例


## 使用snakemake进行样品批量处理

### snakemake的运行逻辑

### 一些可能遇到问题的解决方式

### 批量处理示例


## snakemake & conda


## snakemak & 容器

## 问题解决方案收集

本人在使用snakemake过程中时不时会碰到一些大小问题, 由于这个软件的中文使用者似乎并不多, 英文材料虽然完善但是由于我相关的储备不足, 并不能很快从中找到解决方案, 因此将我碰到的问题&解决方案收集在此, 部分问题的描述可能并不专业, 以后有机会慢慢更正