---
title: Two Bash Tricks
categories: Script
date: 2019-06-05 21:17:03
tags: ["linux", "bash"]
---

While chaining various programs using Bash, you often discover some amazing and useful tricks. Here are two I recently found.

<!-- Abstract part -->
<!-- more -->

The first trick is something I discovered when trying to replace `samtools` with `sambamba`. `sambamba` is a tool designed to replace `samtools` for handling BAM files, and many of its commands are similar to those in `samtools`. However, one subcommand did not support standard output, which meant that you had to generate an intermediate file before reading it for the next step. This was very wasteful of I/O resources. After searching online, I found out that there are three special files in Linux: `/dev/stdin`, `/dev/stdout`, and `/dev/null`. These correspond to standard input, standard output, and discarding, respectively. For programs that do not explicitly support standard input/output, you can try replacing the input/output files with these special files to achieve reading content from standard input or writing output to standard output:

```bash
sambamba merge -t 4 /dev/stdout input1.bam input2.bam \
| sambamba sort -o /dev/stdout /dev/stdin \
| samtools view -Sb -f 30 \
> output.srt.flt.bam
```

`/dev/null` is also very useful. For example, `samblaster` only supports the `-o` option by default, which outputs the main result to standard output. However, I actually need the `-s` or `-u` options and want to process these results further. So, I can use `/dev/null` to discard the main result and then output the needed results to `/dev/stdout` for subsequent operations:

```bash
samtools view -h \
| samblaster \
    --ignoreUnmated \
    -o /dev/null \
    -d /dev/null \
    -s /dev/null \
    -u /dev/stdout \
| sed "s/_[0-9]//g" \
| /soft/bwa-0.7.15/bin/bwa mem \
    -t 30 -k 32 -M -R "@RG\tID:${sam_id}\tLB:${sam_id}\tSM:${sam_id}\tPL:ILLUMINA" \
    /resource/GV/ref_genome/hg38/bwa-0.7.15/hg38.fa /dev/stdin \
| samtools view -Sb -f 256 \
| samtools merge -f output.bam input.bam /dev/stdin
```

Another trick is to mask command results as files for input. For example, if I have multiple two-column files where one column is a data title and the other is a data value, and I want to merge them into one while retaining the first column's titles, I can do:

```bash
paste file1 \
    <(cut -f 2 file2) \
    <(cut -f 2 file3) \
    <(cut -f 2 file4) \
    > pasted_file
```

In the above operation, the data from `file{2..4}` is intercepted using `cut -f 2` and then masked as files to be input into the `paste` command. This avoids using intermediate files and makes the process more efficient.
```