---
title: 'Marp的使用'
categories: Script
date: 2022-06-22 15:40:28
tags: ['Marp', 'VSCode', 'Markdown', '演示文档', 'Python']
---

日常工作中总有需要做一个演示文档展示项目进度或项目结果的时候, 但做PPT其实是个比较花时间的工作(强迫症总是控制不住自己调字体和位置的手), 然后我每次展示的东西都只有一些简单的文字和图表, 从来没用过PPT那些高端大气上档次的元素和功能, 所以找了Marp来偷懒...


<!-- 摘要部分 -->
<!-- more -->


Marp是一个用markdown语法生成演示文档的项目, 用最基本的markdown语法即可快速生成一份还能看(其实有点丑)的掩饰文档出来, 支持生成pptx, pdf以及网页, 用了有一个月了体验还是很好的(看的人可能会觉得眼瞎?)

## Marp_VScode 基本使用
Marp项目下有多个可用程序, 但是用起来最方便的当然还是VScode下的插件了, 不需要特别设置, 安装即用. 


![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_20220622_1655884337.jpg)

安装完成后打开markdown文件时编辑器右上角就会多出一个按钮, 可以切换当前文件为Marp模式

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Set_20220622_1655884577.jpg)


之后就可以开心的写的了... 因为都是基本markdown语法, 除了用`---`将每页的内容分隔开, 也没啥要注意的. 有啥看[官方文档](https://marpit.marp.app/markdown)就行了

## 用例备查

### 插入图片左/右置

我的展示中经常需要放一张图然后对图片做一些说明, 这种事后图片在一侧, 文字在一侧的方式看上去比较合理, 但是markdown的基本语法显然没这么个东西... 可以通过Marp的背景图片语法来实现.

这里的`contain`是让图形自动调整大小, 否则一般都会超出幻灯范围, 例子用了右置, 用`left`也可以丢左边, 然后上下的话似乎`contain`就没效果了...

```markdown
---
marp: true
---

## 右侧是一张统计图
- 图是乱找的

![bg contain right](https://seaborn.pydata.org/_images/grouped_barplot.png)
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide1_20220622_1655885449.jpg)

### 脚标和页码

在开头增加两个设置项就能打开脚标和页码, 不过位置好像用markdown的语法没法调

```markdown
---
marp: true
footer: 20220610

---
## 展示页码用
```

![](https://raw.githubusercontent.com/silenwang/Gallary/master/2022/06/upgit_Marp_VScode_Slide2_20220622_1655885767.jpg)
