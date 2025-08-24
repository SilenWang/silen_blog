---
title: Using Snakemake for Workflow Control in Unconventional Ways
categories: Coding
date: 2018-11-17 12:47:57
tags: ['snakemake']
---

<!-- Summary -->
Snakemake is indeed a very useful tool for workflow development and management. However, in certain scenarios, it can also bring some issues. Coincidentally, I discovered a rather unconventional way of using it: only use Snakemake's dependency handling and task management, while generating scripts separately.
<!-- more -->

The reason was that when writing the workflow, I had to place an actual executable script at a specific location so that users could read these scripts and modify specific steps before resubmitting for execution.

The original intention was good: allowing users who are not familiar with Snakemake to easily modify analysis parameters for specific samples. However, the problem is that Snakemake does not support exporting actual executable scripts as files; it can only output the execution content directly in the log using the `-p` parameter. Moreover, before finding a solution, I was unaware that wildcards could be used in cluster submission commands... So, after several twists and turns, I found this rather unconventional solution.

1. Divide the entire workflow into several parts and write programs to generate all required scripts.
2. When submitting to SGE, use the `-sync y` parameter with `qsub` so that the submission does not end immediately but waits for the completion/failure of the submitted job. If the task fails, it will return a non-zero value.
3. Use Snakemake to build workflow dependencies as usual, specifying `input` and `output` normally. Add the execution script to `input`, and only write the submission command in `shell:`.

This way, we can solve the problem of not being able to directly generate scripts by generating them ourselves and submitting them.
At the same time, we can also take advantage of Snakemake's dependency resolution and task management features for progress monitoring and breakpoint continuation.
However, this approach also digs a huge pit: the script execution content and file names are not specified within Snakemake. If they do not match, the workflow will fail, and there is no good method to locate such issues; it relies entirely on the writer's caution.

For example, you write a script `test.sh` like this:

```bash
#$ -sync y
cat test.txt | cut -f 2 > tset.col2.txt
```

And your Snakemakefile is like this:

```python
rule all:
    input:
        test.col2.txt

rule test:
    input:
        {sample}.txt
    output:
        {sample}.col2.txt
    shell:
        '''
        qsub test.sh
        '''
```

When using Snakemake to run, it will never succeed... because the script generates a result file name that does not match the target file name in Snakemake.
By default, since `rule test` is determined to have failed, other result files used for judgment (if they exist) will be deleted...

In other words, using this method requires you to ensure consistency between the target file names. Excluding this drawback, this approach can fully leverage Snakemake's advantages in task management while avoiding the time required to understand the wildcard mechanism when initially getting started with Snakemake.
