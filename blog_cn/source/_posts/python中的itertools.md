---
title: python中的itertools
categories: Script
date: 2019-02-13 00:37:27
tags: ["Python", "迭代器", "生成器", "itertools", "groupby", "tee", "chain"]
---

原来写脚本的时候很多实用的功能都没有用过, 这次学习了一点点, 做个记录
<!-- 摘要部分 -->
<!-- more -->

本次涉及到对大文件进行处理, 但是大部分时候是单行处理的, 所以想到用itertools以生成生成器. 这样之后也许可以将生成器分割成小部分, 应用到并行运算中.

本次用到的itertools中的三个函数:

1. groupby
2. tee
3. chain

- groupby用于对迭代对象进行聚类, 对应实例是UMI分组时, 需要对比对位置接近的reads进行判断, 将比对情况认为一致的reads放一起进一步判断UMI情况

```python
# 后面是一个以iter_obj的单个对象作为输入的函数, 以函数的返回值作为判断依据
grp_generator = groupby(iter_obj, lambda x: x[0])
# 需要注意的是, groupby只将临近的对象聚在一起, 所以有必要的话自行排序
```

- tee则很简单, 如果迭代对象是生成器的话, 由于生成器只能用一次的, 对于要用到两次及以上的生成器就需要用tee函数做个副本出来

```python
iterable_fork_1, iterable_fork_2 = tee(iterable, 2) # 2可以不写, 默认就是2
```

- chain则是酱多个迭代对象的接过连起来, 我本次具体是用来展平嵌套结果:

```python
flat_iter = chain.from_iterable(ori_iter)
# [[1,2],3] -> [1,2,3] 当然实际生成的是生成器, 要自己list(flat_iter)获取完整结果
```
