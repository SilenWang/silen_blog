---
title: Graphviz的使用
categories: Others
date: 2019-04-25 21:31:53
tags: ['流程图','Graphviz','可视化工具','软件开发']
---

进行软件开发时, 不论是为了展示或者是自己梳理需要都会要用到流程图, 如果自己开个绘图板或者word来画呢也不是不可以, 但是如果图一单复杂了, 或者图是不是会需要更改的话, 就很花时间了, 因此一个类似markdown这样可以让人专注与内容而不是形式的工具就非常重要了.

<!-- more -->

Graphviz是一个脚本化做图的工具, 其可以制作的图形很多, 我目前只使用了流程图. 作流程图需要编写一个`.dot`格式的文件, 这个文件的形式如下:

```Graphviz
digraph workflow{
    // 图形的总定义部分
    label = "Workflow For Mutation Call"
    node[shape = box]

    // 流程图元素的定义部分
    raw [label="Raw Bam"];
    clean [label="UMI Dedup Bam"];
    align [label="Split Bam"];
    vcf [label="Disc Bam(Prototype)"];

    // 流程图关系部分
    raw -> clean [label="fastp"];
    clean -> align [label="bwa"];
    align -> vcf [label="GATK Mutect2"];
}
```

上面的代码我分成了三个部分:
- 第一部分是图形中元素的默认设定, 比如`label`是图形的标题, `node`是流程图的各个节点, 我把默认的节点图形设定为方框(默认是圆)
- 第二部分是图形中节点的设定, 即定义有哪些节点, 然后在后面的`[]`中对节点进行自定义, 比如`label`就可以指定显示的文字
- 第三部分是流程图的顺序关系, 比如`raw -> clean`就从`raw`到`clean`画了一个箭头, 同样后面可以用`[]`自定义

完成上面代码编写后, 就可以调用`graphviz`的`dot`命令生成图片了:

```bash
dot -Tpng workflow.dot > workflow.png
```

- ![完成图形](https://raw.githubusercontent.com/SilenWang/Gallary/master/workflow.png)

vscode下有相应的插件, Linux上安装插件后就可以边写边预览流程图, 非常方便.
