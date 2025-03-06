---
title: Graphviz Usage
categories: Others
date: 2019-04-25 21:31:53
tags: 'graphviz'
---

When developing software, whether for presentation or to organize thoughts, using flowcharts is often necessary. Instead of using a drawing board or Word, which can be time-consuming when the diagram becomes complex and needs frequent updates, a tool like Markdown that allows you to focus on content rather than form is very important.
<!-- Abstract part -->
<!-- more -->

Graphviz is a script-based graph visualization tool that can create many types of graphs. Currently, I only use flowcharts. To create a flowchart, you need to write a `.dot` formatted file, which looks like this:

```Graphviz
digraph workflow{
    // Default definition section for the entire graph
    label = "Workflow For Mutation Call"
    node[shape = box]

    // Definition section for nodes in the graph
    raw [label="Raw Bam"];
    clean [label="UMI Dedup Bam"];
    align [label="Split Bam"];
    vcf [label="Disc Bam(Prototype)"];

    // Relationship section for the flowchart
    raw -> clean [label="fastp"];
    clean -> align [label="bwa"];
    align -> vcf [label="GATK Mutect2"];
}
```

I divided the above code into three parts:
- The first part sets default properties for elements in the graph, such as `label` which is the title of the graph, and `node` which defines the nodes in the flowchart. I set the default node shape to a box (the default is a circle).
- The second part defines the nodes in the graph by specifying what nodes there are and then customizing them in the brackets, such as `label` to specify the displayed text.
- The third part defines the order of relationships in the flowchart, such as `raw -> clean` which draws an arrow from `raw` to `clean`. Similarly, you can customize it with brackets.

After writing the above code, you can generate images by calling the `dot` command of Graphviz:

```bash
dot -Tpng workflow.dot > workflow.png
```

- ![Completed graph](https://raw.githubusercontent.com/SilenWang/Gallary/master/workflow.png)

There is a corresponding plugin for VSCode, and on Linux, after installing the plugin, you can preview flowcharts as you write, which is very convenient.
```