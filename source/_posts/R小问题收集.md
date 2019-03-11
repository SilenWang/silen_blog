---
title: R小技巧收集
categories: Script
date: 2018-09-07 12:01:44
tags: ['R']
---

用R总是会碰到各种问题, 收集一下
<!-- more -->

# 数据框操作

## 取单列/行不返回向量

- 使用`df[]`的方式对数据框取子集时, 若取出的是单列/行, 或者删除后只剩下但列/行的话返回的对象会是向量而不再是数据框, 需要的时候要加参数指定以返回数据框
- 本方法对matrix一样有效

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

# 其他
## R中的Python Dict替代

- 许多网上的资料会使用list配合match函数实现类dict功能, 但其实R中的向量是可以被命名的, 通过命名向量即可实现key-value的映射, 十分方便. 不过这么做的效率是否好就不得而知了..

```r
a <- c(1, 2, 3)
names(a) <- c("a", "b", "c")
tmp <- a["a"]
print(tmp)
```