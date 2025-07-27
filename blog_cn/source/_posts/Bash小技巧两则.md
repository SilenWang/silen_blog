---
title: Bash小技巧两则
categories: Script
date: 2019-06-05 21:17:03
tags: ["Linux", "Bash", "命令行", "输入输出重定向", "进程替换"]
---

使用Bash串联各个程序的过程中总是能发现一些神奇好用的技巧, 这里记录两个.

<!-- 摘要部分 -->
<!-- more -->
第一则是我前些时候尝试使用sambamba取代samtools的时候发现的, sambamba是用来取代samtools处理bam文件的一个工具, 其很多命令都跟samtools非常接近, 但是当时有一个子命令却不支持标准输出, 导致要先生成中间文件然后再读取中间文件进行下一步, 非常浪费IO资源, 所以我搜索了一下有没啥好的解决方案. 不搜不知道, 一搜才知道, 原来Linux下面是有三个特殊文件的, 分别是`/dev/stdin`, `/dev/stdout`, `/dev/null`, 分别对应标准输入, 标准输出和丢弃. 对于没有明确说明支持标准输入/输出的程序, 可以尝试将输入/输出文件替换为上述特殊文件, 这样就能实现从标准输入读取内容或者将接过输出到标准输出了:

```bash
sambamba merge -t 4 /dev/stdout input1.bam input2.bam \
| sambamba sort -o /dev/stdout /dev/stdin \
| samtools view -Sb -f 30 \
> output.srt.flt.bam
```

`/dev/null`同样很有用, 比如`samblaster`, 默认只支持`-o`, 也就是主接过输出到标准输出, 但是我其实需要`-s`或者`-u`的接过, 并且我需要对这些接过做另外的处理, 那么就可以使用`/dev/null`将主要结果丢弃, 然后将需要的结果输出到`/dev/stdout`然后继续后续操作即可:

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

另外一个技巧是将命令结果伪装成文件进行输入. 比如我现在有多个两列的文件, 一列是数据标题, 一列是数据值. 我现在需要把多个文件合并成一个, 同时还想保留第一列的标题, 那么就可以:

```bash
paste file1 \
    <(cut -f 2 file2) \
    <(cut -f 2 file3) \
    <(cut -f 2 file4) \
    > pasted_file
```

在上面的操作中, `file{2..4}`的数据被`cut -f 2`截取出来, 然后伪装成文件, 输入给了`paste`命令, 这样避免了使用中间文件, 更加快捷高效.
