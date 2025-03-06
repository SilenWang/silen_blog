---
title: SGE System Usage Notes
categories: Bioinformatics
date: 2018-08-30 09:36:38
tags: ['sge']
---

When using SGE at a new company, the approach is somewhat different. Here are some notes.

<!-- more -->

# Script Writing

After writing the script, it is not specified when submitting; instead, it is written into the script to be submitted and then directly `qsub shell.sh`.

```bash
#！/bin/bash
#$ -S /bin/bash          // Indicates that this script is for bash
#$ -V                    // Passes all environment variables from the current command
#$ -cwd                  // Sets the current path as the working directory
#$ -N WorkName           // Task name
#$ -o WorkName.log       // Task output log file name
#$ -j y                  // Specifies whether the standard error stream of the task is merged with the standard output stream [yes] no

shell script
```

Submitting a script containing the above content is essentially equivalent to:

```bash
qsub -S /bin/bash \
    -V \
    -cwd \
    -N WorkName \
    -o WorkName.log \
    -j y \
    shell script
```

# qlogin Usage

For testing purposes, you can use `qlogin` to log in to a compute node. Then, the operations performed on the compute node environment are consistent with submitting a script using `qsub`.
```