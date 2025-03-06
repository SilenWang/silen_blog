---
title: Using GNU Parallel in Bash for Parallel Execution
categories: Script
date: 2019-03-31 23:07:26
tags: ["bash", "parallel computing"]
---

I recently discovered GNU Parallel, a simple and useful parallel program implemented in Perl, which can easily execute various commands in Bash in parallel. This is particularly useful for tasks that can be divided into independent subtasks.

<!-- Abstract part -->
<!-- more -->

GNU Parallel allows you to replace loops with parallel execution by specifying the command to be executed in parallel. For example, consider a script that filters VCF files based on chromosomes:

```bash
for chr in {{1..22},X,Y};do
    bcftools view -f PASS -r chr${chr} sample.vcf > sample.chr${chr}.flt.vcf
done
bcftools merge sample.chr{{1..22},X,Y}.vcf > sample.flt.vcf
```

Using GNU Parallel, you can rewrite it as follows:

```bash
echo {{1..22},X,Y} | tr " " "\n" | parallel 'bcftools view -f PASS -r chr{} sample.vcf > sample.chr{}.flt.vcf'
bcftools merge sample.chr{{1..22},X,Y}.vcf > sample.flt.vcf
```

In the first example, the script loops through each chromosome and filters the VCF file. In the second example, the filtering is done in parallel, followed by merging. Of course, this example doesn't make sense because both methods essentially do the same thing.

By default, GNU Parallel uses all available cores, so you can specify a maximum number of cores using `-j` or set a load threshold with `--load 80%`.

With this tool, any command that supports positional parameter analysis for regions can be quickly converted to parallel execution, significantly speeding up the analysis process. This makes writing Snakemake workflows more convenient...

I will continue to explore other uses of this program in future posts...
