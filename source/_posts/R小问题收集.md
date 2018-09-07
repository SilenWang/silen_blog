---
title: R小问题收集
categories: Script
date: 2018-09-07 12:01:44
tags: ['R']
---

用R总是会碰到各种问题, 收集一下
<!-- more -->

# 数据框操作

## 取单列/行不返回向量

- 使用`df[]`的方式对数据框取子集时, 若取出的是单列/行, 或者删除后只剩下但列/行的话返回的对象会是向量而不再是数据框, 需要的时候要加参数指定以返回数据框

```r
df[c(1, 2, 3), -2, drop = F]
```

## 合并数据时以row.names为key

- 大多中文材料对合并时的描述都是采用某一列数据作为key进行列合并, 其实也可以以row.names作为key, 但是合并后row.name就变成一列数据了, 需要再次指定row.names

```r
# 分两步进行
merge_data <- merge(x, y, all.x = T, by = 0)
merge_data <- data.frame(merge_data), row.names = 1)
# 合并成一步
merge_data <- data.frame(merge(x, y, all.x = T, by = 0), row.names = 1)
```