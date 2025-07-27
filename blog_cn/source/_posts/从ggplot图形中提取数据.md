---
title: 从ggplot的箱线图中提取离群值并进行去除
categories: Script
date: 2020-01-05 17:56:20
tags: ['R语言', 'ggplot2', '数据可视化', '箱线图', '离群值处理', 'R', 'ggplot']
---

今天碰到一个做图问题, 学到了一个新的小技巧, 记录一下
<!-- 摘要部分 -->
<!-- more -->

ggplot提供了很多常用的图形的绘制函数, 相当方便, 且图形可以按照各种各样的分组进行绘制, 但是有的时候我有些特殊要求, ggplot无法满足, 所以需要自行进行一些操作.

比如我今天需要绘制一个分组的箱线图, 每一组的数据中都有自己的离群值, 这些离群值会使图形挤在某一测, 难以观察, 因此需要去除. `geom_boxplot`虽然有离群值处理的参数, 但是很可惜是这些参数并不会把离群点从数值中去除, 而是不在图形上显示, 这么一来该挤一起的图形还是挤一起, 达不到效果.

网上能查到的关于离群值去除的方案最多是利用内置`boxplot`函数里面`$out`的内容提取离群数值并进行去除, 可惜不巧的是这这种方案似乎不能分组给出离群值结果, 它给的是总体离群值, 如果需要分组去除的话还需要自己手动拆分数据集然后挨个处理.

那么有没有更偷懒一点的方式呢?

还真被我发现了...

ggplot的图形和内置函数一样, 是可以提取实际绘图用数据的:

```r
plot <- ggplot(data=data,mapping=aes(x=data$label, y=data[,tag])) +
    geom_boxplot()
plot_data <- layer_data(plot)
```

`layer_data`函数是提取的方式之一, 使用`build_plot(plot)$data`效果一样. 

`plot_data`中就存储这ggplot绘制图形时使用的一切数据. 其中我们需要的部分存储在`outliers`中. 并且顺序是按照图形上展示的顺序排列的. 因此只要在前期绘图时对分组手动指定factor level, 就可以对应从图形中拿到特定组的离群值.

然后利用ggplot绘图时自动去`NA`的特性, 将对应列的离群值以`NA`替换, 就能轻松的创造一个去过离群值的数据集, 然后再次画图就能解决问题啦~

完整的绘图代码与图形对比:

```r
library('ggplot2')
data <- read.csv("plot_data.tsv", sep="\t", stringsAsFactors=F, na.strings = '.')
# 对分组数据自定义等级, 后续对应必须
data$grp <- factor(data$grp, levels=c("Grp1","Grp2"))
# 准备一个vector作映射
classV <- c(1,2)
names(classV) <- c("Grp1","Grp2")
# 第一轮, 利用ggplot获取离群值
plot <- ggplot(data=data,mapping=aes(x=data$grp, y=data[,'value'])) +
    geom_boxplot()
print(layer_data(plot))
outlier_data <- layer_data(plot)['outliers']
# 利用获取数值改变原数据集的离群值为NA
pdata<-data
for (grp in names(classV)) {
    pdata[pdata$grp==grp,  'value'] <- replace(pdata[pdata$grp==grp, 'value'], pdata[pdata$grp==grp,  'value'] %in% outlier_data[classV[grp], ][[1]], NA)
}
# 准备比较组向量, 利用变更后的和数据作图
comb_list <- list()
comb <- combn(c("Negative", "Unkown","Weak", "Medium","Strong"),2)
for (id in seq(dim(comb)[2])) {
    comb_list[[id]] <- comb[, id]
}
pplot <- ggplot(data=pdata,mapping=aes(x=pdata$grp, y=pdata[,'value'])) +
geom_boxplot()
```

![display_boxplot.png](https://raw.githubusercontent.com/SilenWang/Gallary/master/display_boxplot.png)
