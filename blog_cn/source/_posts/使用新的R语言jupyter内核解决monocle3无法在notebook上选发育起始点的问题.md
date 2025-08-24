---
title: 使用新的R语言jupyter内核解决monocle3无法在notebook上选发育起始点的问题
tags:
  - 单细胞分析
  - 拟时间分析
  - monocle3
  - jupyter
  - xeus-r
categories: Bioinformatics
date: 2024-06-23 01:47:03
---


看一看这个草稿的创建时间，居然是2024年06月了... 现在都2025年6月了，突然理解了为什么那么多UP主会成为鸽子，挖坑一时爽，填坑火葬场啊... 这篇是我位数不多的纯生物信息内容了...

问题起因是去年要运行 monocle3 进行拟时间分析，但是分析到最后发现一个很蛋疼的问题，monocle3 中的拟时间轨迹起点，是需要分析者来手动指定的，R代码执行过程中会自动开启浏览器，用户需要在浏览器打开的临时网站中指定起点，然后关闭页面，分析继续进行。

但是jupyter的`irkernel`内核并不支持这一特性。也就是说，我无法在juputer notebook中直接完成分析。这一问题在[2019年就被提出过](https://github.com/cole-trapnell-lab/monocle3/issues/179)，但是直到我需要进行分析的2024年，5年过去，依然没有解决方案...

<!-- more -->

不过我当时还真偶然发现了解决方案，因为github给我推荐了一个项目，[jupyter-xeus](https://github.com/jupyter-xeus/)，这个项目旨在重新开发一系列Jupyter内核，其中就包括新的R内核，[xeus-r](https://github.com/jupyter-xeus/xeus-r)，这个新的内核支持交互模式，可以在运行 monocle3 的时候生成一个链接，用户点进链接完成选点后，再退出页面就能继续分析了。

`xeus-r` 在conda上就有，用 conda / mamba / pixi 就能快速安装。安装完后，在jupyter的界面会看到对应的`xr`内核，选择使用就好

![xeus-r](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/06/upgit_20250615_1749994103.png)

需要执行的代码如下：

```r
library(monocle3)
cds <- readRDS("/Your/rds/contain/monocle3/object.RDS")
options(browser="firefox")
cds <- order_cells(cds)
```

运行后在单元格下面就可以看见需要打开的链接，打开后，选择起点，再后关闭即可。

![monocle3](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/06/upgit_20250615_1749997047.png)
