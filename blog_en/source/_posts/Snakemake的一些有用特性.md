---
title: Some Useful Features of Snakemake
categories: Script
date: 2019-03-06 23:41:40
tags: ['Python', 'snakemake']
---

After using snakemake for some time, I found that it has many practical features to facilitate daily analysis. Here is a record and整理.

<!-- more -->

## Resource Limits

In the `rule`, you can use the `resources` keyword to limit resources. This is very useful when you need to use a GPU (only one GPU) or have limited disk performance (cannot perform too many high IO operations at the same time). For example, the following setting:

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

In the above workflow, the sort process uses 8 threads simultaneously and consumes a lot of IO. Therefore, its IO usage is set to 8, while the call operation sets `IO=1`. When running, you can control the number of high IO tasks by setting the maximum IO value: `snakemake --resources IO=16`.

## Task Groups

When using a cluster, in certain situations, grouping related tasks to be submitted to the same node can save on IO (prevent来回 transmission of intermediate files between nodes, roughly speaking...), and the keyword used is `group`. Tasks within the same group will be submitted to the same node for execution at once (the node will start another snakemake instance to control them).

The actual usage is very simple, so I won't write an example. However, there is one thing to note: after tasks are submitted to a node using this feature to utilize the local disk of the node, although there is a snakemake controlling them, when the task submission ends, the local snakemake process will once again check if the input and output meet the requirements. This point needs to be noted. For example, I used this function to take advantage of the local disk of the compute node, but these disks are limited to the compute nodes and cannot be accessed from the login node. This will cause the login node's process to fail when checking if the tasks are completed because some intermediate results cannot be found.
```