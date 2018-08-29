---
title: Pandas操作笔记
categories: Script
date: 2018-08-29 12:16:24
tags: ['python', 'pandas']
---

最近用pandas比较多, 把常用的东西摘抄一下备查

<!-- more -->


# 数据读取

- 基本读取直接read_table就可以了, 主要是涉及到一些参数

```python
import pandas as pd
df = pd.read_table(file_path, sep="\t", index_col=0, header=0)
# index_col等同R里面的row.names, 都是把特定列作为行标签使用, 并且将该列从数据中去除, 如果不指定则会生成0-length的数字作为标签
# header函数与R里的逻辑不太一样, 默认是header=0, 即将读取的第一行作为表头, 如果不要表头的话用header=None, 如果制定别的行为表头, 则表头行以上的数据会被丢弃
```
- 在处理很大的数据时, 为防止爆内存, 需要生成迭代器分块读取文件

```python
import pandas as pd
# reader是以一个生成器
reader = pd.read_table(file_path, sep="\t", iterator=True, chunksize=1000))
# iterator / chunksize这两个参数指定一个就会生成迭代器, 其实如果指定了chunksize可以不写iterator了

# 如果生成reader时指定过chunksize, 那可可以直接调用get_chunk()
# 获得特定数量的行, 否则就要写get_chunk(num), 这个我还没详尽测试过
# 注意, 这里的行数是数据的行数, 不包括表头, 因为生成reader的时候
# 默认以第一行作为表头了, 所以后面生成的df都是带表头的
chunk = reader.get_chunk()
df = pd.DataFrame(chunk)
```

# 子集选取

- 与R类似, pandas也支持多种子集选取方式, 部分选取逻辑也跟R非常类似, 但是有些东西与R并不一样, 一般来说在pandas里能实现的方式R都可以, 但是反过来则不然, 这里只记录与R中选取逻辑类似的方法

```python
import pandas as pd
df = pd.read_table(file_path, sep="\t")
# 依据标签选取行列使用.loc[row, col]
df = df.loc[['a', 'b', 'c'], ['d', 'e', 'f']]
# 注意与R不一样的是, 单取列时, 前面不可以留空, 单取行则可以后面留空
df = df.loc[:, ['d', 'e', 'f']]
df = df.loc[['a', 'b', 'c'], ]
# 与.loc对应有一个.iloc, 这里使用的行/列序号, 而不是标签, 注意需要从0起算

```

# 行列删减

# 遍历处理