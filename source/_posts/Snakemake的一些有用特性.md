---
title: Snakemake的一些有用特性
categories: Script
date: 2019-03-06 23:41:40
tags: ['python', 'snakemake']
---

在使用snakemake一段时间后, 发现其中确实有很多实用的特性来方便平常的分析, 做个记录整理

<!-- more -->

## 资源限制

在`rule`内可使用`resources`关键字来限定资源, 这在需要用到GPU(只有一个GPU)或者磁盘性能有限(不能同时进行太多高IO操作)的时候非常有用. 如下面的设定

```python
rule all:
    expand("{sam}.vcf", sam=["sam1", "sam2", "sam3", "sam4"])

rule samtools_sort:
    input: 
        rules.samtools_flt.output
    output:
        "{sam}.flt.srt.bam"
    resources:
        IO=8
    shell:
        "samtools sort -@ 8 {input} > {output}"

rule samtools_call:
    input: 
        "{sam}.bam"
    output:
        "{sam}.vcf"
    resources:
        IO=1
    shell:
        '''
        bcftools mpileup -Ou -f reference.fa {input} \\
            | bcftools call -mv -O -o {output}
        '''
```

上面的流程中, sort过程8线程同时进行, 相当消耗IO, 因此设定其使用的IO为8, 而call操作设`IO=1`, 之后在运行的时候设定最大IO数值就可以达到控制高IO任务数目的目的: `snakemake --resources IO=16`

## 任务组

需要使用集群时, 在特定的情况下, 将相互关联的任务投递到同一个节点会更节省IO(防止节点间来回传输中间文件, 大概吧...), 使用关键字是`group`, 同一个组内的任务会被一次性的投递到同一个节点中去执行(节点上会再启动一个snakemake用于控制).

实际使用很简单所以就不写例子. 但是需要注意一个问题, 同组同任务投递上节点后, 虽然有snakemake进行控制, 但是在投递任务结束后, 本地的snakemake进程是会再一次检查input和output是否满足的, 这一点要注意. 比如我使用这个功能是为了利用节点的本地磁盘, 但是这些磁盘是计算节点限定的, 登录节点并不能访问. 这就会导致登录节点的进程检查任务是否完成时, 部分中间结果找不到而导致任务执行失败.