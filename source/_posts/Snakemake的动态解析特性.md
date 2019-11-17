---
title: Snakemake的动态解析特性
categories: Script
date: 2019-11-17 23:52:31
tags: ['Python', 'snakemake']
---

特么我终于成功使用这个特性了...我真是太南了...

<!-- 摘要部分 -->
<!-- more -->

在实际的流程中, 我们的流程可能会有一些"动态"的要求. 比如对fq进行比对, 或者对bam进行变异检测, 对应的比对/检测软件可能并没有完善的多线程机制, 同时原始数据文件过分的大, 直接进行运算非常耗时. 这时候, 使用`源文件拆分->多进程并行->结果文件合并`就是一种非常实用的策略. 那么, 怎么拆呢? 我可以将文件按固定的方式来拆, 比如变异检测时经常会将bam按照染色体拆分开. 但是这种固定的拆分模式难免会出现效率不尽如人意的问题. 比如不同染色体的数据量是完全不一样的, 按染色体并行, 速度必然受数据量最大的那个染色体的限制. 但是如果不按照固定策略拆, 生成的拆分文件个数很可能是不固定的, 这与`snakemake`先解析要做什么再分配运行的逻辑冲突, 因此开发者在较新的版本中引入了动态解析机制(`Data-dependent conditional execution`).

这个机制的核心是`checkpoint`, 检查点是一种特殊的`rule`, 其在运行时与`rule`没有太大不同, 但是, 在流程运行到检查点后, 整个依赖图(`DAG`)会重新进行解析. 设计者的用意很明显: 将产生不定结果的规则设定为检查点, 在不定的文件实际产生后重新决定需要进行的步骤.

由于这个动态机制十分特殊, 因此检查点使用起来也比较麻烦, 无法像`rule`一样简单的使用. 

首先`checkpoint`必须配合`function input`使用, 因为虽然加了新机制, 但是这个程序的运行逻辑依然是决定做什么然后分配. 而决定做什么这一点不是程序能自己决定的, 需要我们手动设计, 因此必须使用一个函数去获取检查点生成的文件, 然后将这些文件交给下游`rule`去解析. 

我这有一个比官网那个稍微好理解点的例子:

```python
from os import path
from glob import glob

rule all:
    input:
        "{sample}.all.stat".format(sample="NL180614-2_S12_L001_R1_001")

checkpoint split:
    input:
        "{sample}.fastq.gz"
    output:
        touch("{sample}.splt.done")
    shell:
        "split -l 1000000 {input} {wildcards.sample}.split. "


rule wc:
    input:
        "{sample}.split.{split}"
    output:
        "{sample}.split.{split}.stat"
    shell:
        "wc -l {input} > {output}"


def get_split_files(wildcards):
    split_output = checkpoints.split.get(**wildcards).output[0]
    split_dir = path.split(split_output)[0]
    split_files = glob(f"*.split.*")
    tar_files = [f"{f}.stat" for f in split_files]
    return tar_files


rule cat:
    input:
        get_split_files
    output:
        "{sample}.all.stat"
    shell:
        '''
        cat {input} > {output}
        '''
```

这个例子实现的是使用`split`将fq文件拆分成多个子文件, 分别计数行数然后将结果合并成一个文件. `checkpoint`的部分是拆分, 同时这个`rule`产生了一个不交由下游解析的结果, 也就是完全依赖无关. 真正有用的其实是`checkpoint`本身. `get_split_files()`函数中使用了`checkpoints`对象, 并且使用了其下`split`这个检查点的内容, 因此以这个函数作为输入的`cat`规则上游依赖是`checkpoint split`, 同时其也不会进行真正的依赖解析, 而是等待它的`checkpoint split`完成后再解析(也就是动态解析).

官方文档的例子基本也是这个原理, 只不过其`checkpoint`生成了目录结果, 以方便重解析生成的动态文件.

明天有空的话写一个实际的拆分后比对的例子再更新上来.

以上~