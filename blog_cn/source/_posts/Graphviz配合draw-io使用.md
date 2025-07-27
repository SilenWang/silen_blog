---
title: Graphviz配合draw.io使用
categories: Script
date: 2021-10-20 01:48:35
tags: ['流程图', '图表工具', 'graphviz', 'drawio']
---

自从会了Graphviz后, 经常会用它来做流程图或者简单的示意图, 大部分使用它都是比较好用的, 因为流程图比较小, 同时也只用来示意, 并不用于展示. 但是当流程图较复杂, 同时有美化和展示的需要, 就会比较麻烦了.
<!-- 摘要部分 -->
<!-- more -->

因此我找到了一个折中的办法, 先用graphviz确定草图, 然后使用[graphviz2drawio](https://github.com/hbmartin/graphviz2drawio)来对写好的graphviz进行格式转换, 转换好的草图再导入[drawio](https://github.com/jgraph/drawio)来手动进行编排

实际用下来这个方案也并不是特别好用... 导入的流程图大致框架是有, 但是原来graphviz中的图形不会保留(圆啊方块什么的), 同时部分graphviz里的对象在drawio中并不存在, 所以用graphviz做草图的话, 最多只能保留基本的骨架, 需要继续调整的东西还是挺多的. 所以之后还是要么再学一个使用drawio的语言, 要么找一个别的类似drawio的工具来导入graphviz后进行进一步的编辑.
