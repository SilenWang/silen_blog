---
title: MultiQC简单使用
categories: Bioinfomatic
date: 2019-07-27 21:02:16
tags: 
- NGS
- 质控工具
- 生物信息学
- 数据分析
- MultiQC
- fastqc
---

MultiQC是一个NGS数据的质控工具, 不过跟很多其他工具不同的一点是, 它本身并直接进行数据获取和指标计算, 而是读入各种常见质控工具的结果文件然后做综合展示.
<!-- 摘要部分 -->
<!-- more -->

目前MultiQC支持的工具已经达到了75种, 且支持自行添加更多信息(虽然我还没试过), 支持工具覆盖了碱基读取到变异检测的整个过程, 出具的报告也相当美观, 非常适合用来做最后的数据质量汇总.

MultiQC的使用很简单, 首先要将其他工具产生的结果文件放到一个统一的路径下, 然后运行`MultiQC /path/to/file`, MultiQC就会检索文件夹下的各种质控结果, 然后生成最后的汇总报告了, 比如我这里有一份使用MultiQC读取fastqc结果生成的fastq质控报告:

![](https://raw.githubusercontent.com/SilenWang/Gallary/master/multiqc_report.png)

可以看见整个报告还是很漂亮的, 在报告顶部会有各个指标的汇总表格, 左侧是各个软件下计算指标的导航栏, 根据不同的软件, 具体的内容可能是展示图或者展示表. 在顶部的汇总表格中, 默认只显示部分指标, 如果想更改选择的指标可以点`configure columns`来对需要显示的指标进行调整.

![](https://raw.githubusercontent.com/SilenWang/Gallary/master/multiqc_report_config.png)

然后还有一个比较方便的功能是, 可以点击`plot`按钮自行选择汇总表中的两个指标做散点图, 这个功能在后续查看一些异常的时候可能会特别有用.

![](https://raw.githubusercontent.com/SilenWang/Gallary/master/multiqc_report_plot.png)

简单的介绍就这么多, 更多高级的用法如自定义新的指标和变更样品名称匹配规则以让同一个样品的结果汇聚到一起等等, 之后有了新的成果之后再继续水~
