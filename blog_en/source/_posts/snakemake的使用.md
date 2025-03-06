---
title: snakemake's Usage
categories: Bioinformatics
date: 2018-10-30 12:40:03
tags: ['Python', 'snakemake']
---

Writing workflows is a common task in bioinformatics analysis, and a mature and well-designed tool can greatly improve the efficiency of work.
<!-- more -->

The simplest way to write a workflow is to write all the commands you need into a shell script and specify the order of tasks using `&&` or other methods. Then, you can run it in the background.

Although this approach is easy, its repetitiveness and error resistance are really worrying. Therefore, various workflow management tools have emerged: I used `sjm` at my previous unit, 10x Genomics developed `Martian`, Broad Institute's `WDL`, and I am currently using `snakemake`.

## Introduction to snakemake
`snakemake` is a bioinformatics analysis workflow management tool developed by the Köster lab at the University of Duisburg-Essen in Germany. Its entire code is written in Python, so it can be easily integrated with Python programs (although it also has a complete set of command-line programs).

The features of `snakemake` include but are not limited to:

- Simple and easy to use
- Can directly use Python code
- Built-in distributed computing cluster support (SGE, LSF, etc.)
- Built-in Conda environment deployment support
- Built-in container support (note that this container is not Docker)
- Automatic handling of dependencies and breakpoint restarts
- Built-in testing functionality (benchmarking)
- Same execution logic as `make`, so people with a background in `make` can quickly get up to speed

## Simple Example of snakemake

In a `snakemake` file, the most basic unit is a `rule`, which represents a step in the workflow. Under a `rule`, you can define input files, output files, and the commands to be executed. Since `snakemake` is based on Python, its syntax format is consistent with Python, using indentation to define control scopes, such as writing a simple `rule`:

```python
rule samtools_sort:
    input: 
        "test.bam"
    output:
        "test.srt.bam"
    shell:
        "samtools sort {input} > {output}"
```

This `rule` named "samtools_sort" executes the step of sorting the input `test.bam` using the `sort` command from `samtools`, and then specifies the output result as "test.srt.bam".

The final `shell:` indicates that this is a shell command, with the command given as a string. Snakemake also supports directly executing Python code; in this case, the keyword needs to be changed to `run:`. Additionally, it supports calling R code and other external scripts. For more details, refer to the [official manual](https://snakemake.readthedocs.io/en/stable/).

## Using snakemake for Batch Processing of Samples

After understanding the basic writing of a `rule`, you can write a series of `rules` and then link them into a workflow. Although `snakemake` is simple and easy to learn, writing a complete workflow that can handle multiple samples still requires some knowledge.

### Special String Formatting in snakemake

From the previous simple example, it can be seen that commands, input files, and output files in the `snakefile` are given as strings. These strings can use similar formatting string methods to reduce hardcoding and increase applicability.

For example, in the `"samtools sort {input} > {output}"` from the previous example, `{input}` will be replaced with **the** string defined after `input` within this `rule` (if there are multiple inputs, they will be listed separated by spaces). Similarly, we can call any content defined within a `rule` to fill in text templates and generate the final commands or strings that need to be executed.

However, one thing to note is that for those who are used to using `format()`, it might look like this:

```python
path="path/to/input"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {path}/{input} > {output}".format(path=path)
```

This is mixing Python's string formatting with snakemake's string formatting. ~~This kind of writing cannot be correctly handled by snakemake; a string can either use the built-in formatting method or the original formatting method of Python. The above example, using `format()`, will be processed according to Python's original method, filling in the string based on the information given in `format()`. Unspecified `input` and `output` will be left empty.~~ Snakemake can handle this situation; however, the writing style needs to be changed:

```python
path="path/to/input"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {path}/{{input}} > {{output}}".format(path=path)
```

For those who are more familiar with Python, you might know that `"\{\{\}\}"` represents not escaping the `{}` within a string but treating them as general symbols. Combining snakemake's execution logic, it is known that snakemake will handle the pure string input itself after Python part execution is complete.

### Execution Logic of snakemake Workflows

Snakemake does not have a dedicated syntax for specifying order or dependencies; workflow connections are automatically completed based on input and output file dependencies. Since snakemake defaults to parsing and completing the first `rule`, it is recommended by the official documentation to create a `rule` named "all". In this `rule`'s `input` (note: not `output`), specify all final files. Then, the program will find which `rule` can produce the required file based on the `output` in the snakefile. It continues to parse the dependencies of this rule until it finds the source or determines that there is no source and reports an error, for example:

```python
rule all:
    "test.flt.srt.bam"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {input} > {output}"

rule samtools_sort:
    input: 
        rules.samtools_flt.output
    output:
        "test.flt.srt.bam"
    shell:
        "samtools sort {input} > {output}"
```

In the above snakefile, the workflow execution logic is as follows:

1. The final target is `test.flt.srt.bam`, search for a rule that can produce this file.
2. Discover that `samtools_sort` can produce the required file, and this rule requires input from `samtools_flt`. Therefore, `samtools_flt` executes first.
3. `samtools_flt` needs input `test.bam`, check if the file exists. If it exists, follow the order `samtools_flt > samtools_sort > all` to execute the workflow. If not found, report a dependency error and stop the workflow.

Due to this backtracking execution logic of workflows, it is very easy to implement breakpoint restarts in workflows: simply backtrack through each step's dependencies to check if they are satisfied. If satisfied, start executing from the current step; if not, continue backtracking until finding the source or failing to find one and reporting an error.

Additionally, snakemake checks output file updates and integrity during workflow execution to determine whether rules run successfully. If a rule fails, the result files of that rule are deleted by default, so unless the snakemake process is terminated abnormally, it is unlikely for a step to be incomplete while the workflow continues running.

### Built-in Useful Functions in snakemake

Since snakefiles can directly write Python code, people familiar with Python can easily use list comprehensions and other similar things to quickly generate information needed for batch processing. For those without much Python experience, snakemake also provides some built-in functions to quickly achieve similar functionality. For example:

- expand:
    + A function used to generate string lists, equivalent to combining list comprehensions and formatted strings. By default, it combines multiple replacement variables in pairs, but it can be set to a mode similar to `zip()`, where it sequentially pairs them, which can meet different needs.
    + `samples = expand("{sample_id}.{fq}.fasta", sample_id=['sam1', 'sam2'], fq=['fq1', 'fq2'])`

### Batch Processing of Multiple Samples/Files

Batch processing involves a built-in object in snakemake called `wildcards`. Although I am not very clear on how to translate this term... In general, this object is used for pattern matching with capture effects (regular expressions). For example, an example given in the official documentation:

```python
rule complex_conversion:
    input:
        "{dataset}/inputfile"
    output:
        "{dataset}/file.{group}.txt"
    shell:
        "somecommand --group {wildcards.group} < {input} > {output}"
```

Here, the `input` and `output` of the rule contain `{dataset}` and `{group}`, which are two `wildcard` objects. Each `wildcard` can match almost any character (corresponding to a regular expression of `.+`). If a rule needs a file named `file_path/file.001.txt`, snakemake will find that this file's naming matches the pattern specified in the `output` of `complex_conversion`. Therefore, snakemake considers `file_path/file.001.txt` to be produced by this rule and sets `dataset` to "file_path", `group` to "001", then passes these variables up to `input`, continuing dependency resolution.

Thus, through dependencies, when the requested result files are named according to a certain rule (usually containing sample/project IDs), you can use `wildcard` to match which steps each sample/project needs to perform and automatically generate the complete workflow route before starting execution.

## snakemake & conda

To ensure the repeatability of completed workflows, snakemake includes built-in functionality for using Conda to manage workflow dependencies.
However, this feature is different from what I initially understood.

## snakemake & Containers

In addition to Conda, snakemake supports using containers to manage dependencies. However, this container is not Docker; it is Singularity. I only discovered that Docker is just one implementation of container counting when using snakemake. The key features of snakemake's Singularity compared to Docker are:

- No need for root permissions
- Lower resource consumption

This project started migrating from C++ to Go in version 3.0, and although I overcame network pressure to install it, due to distribution issues, I cannot use it normally at present, so I do not have practical experience with it.

## Problem Solutions Collection

Throughout my use of snakemake, I occasionally encounter some large problems. Since the number of Chinese users for this software is not many, although the English materials are comprehensive, due to my lack of relevant knowledge, I cannot quickly find solutions from them. Therefore, I have collected the problems and their solutions here. Some descriptions of the problems may not be professional; in the future, I will gradually correct them.
