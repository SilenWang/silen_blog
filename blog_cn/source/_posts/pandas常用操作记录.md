---
title: pandas常用操作记录
categories: Coding
date: 2020-07-26 00:36:49
tags: ['数据处理', '缺失值处理', 'Pandas', 'Python']
---

<!-- 摘要部分 -->
<!-- more -->
因为现在的工作有大量数据统计/处理要求, 且论错误处理及第三方扩展性, Python方便得多, 因此目前我数据处理也开始从R移动到Pandas了. 在此记录一些实际工作中常用的操作及需要注意的问题:

## 1. apply相关

### 1.1. 基本使用

同Pandas与R一样有`apply`函数, 不过仅此一个, 没有`apply`家族. `apply`本质上就是对表格进行按行循环, 而之所以用`apply`而不是直接`itterrow`, 则是出于并行的考虑. 有非常多的第三方模块可以直接将apply操作并行化, 因此最初写代码时就按照`apply`的写法来, 可以省下日后并行化时重写函数的过程.

当然另外一点就是将主处理逻辑和细部的详细执行内容分开, 一定程度上可以提高代码可读性, 不过个人感觉还是个习惯问题... 如果读代码的人从不这么写... 那用`apply`其实是提高读代码人的阅读难度...

目前我使用的方式一般是实际计算的部分写成通用函数, 然后套`apply`时使用匿名函数传参数.

```python
def plus(x, y):
    return x+y

data['applyed'] = data['2apply'].apply(lambda x: plus(x['x'], x['y']))
```

当然, 如果同时要对同一行做很多行内操作, 那就单独写个`func_apply()`的函数, 防止反复对数据遍历.

### 1.2. 返回多列结果

`apply`的常规用法是返回一个数值形以形成一列数据, 但是实际应用中我其实会有返回多列的需要, 因此从网上找了一种解决方案:

```python
def plus(x, y):
    return x+y, x-y


data['x+y'], data['x-y'] = zip(*data['2apply'].apply(lambda x: plus(x['x'], x['y'])))
```

## 2. Groupby相关

### 2.1. 多行并为单行数据

这个操作主要是用来将多行的数据统计为一行结果, 我印象中R里面有专门的处理方式, pandas应该也有...然而我一时没找到, 所以目前自己想了一种处理, 这种方式本质上是按待处理单位获取数据子集, 根据子集统计后生成新的数据行填入新表. 这种方式目前来看效率堪忧, 且因为本质是用了for循环, 并行化就没有`apply`那么简单了.

```python
result_list = []

for k, sub_data in data.groupby(by=['Grp1', 'Grp2']):
    grp1, grp2 = k
    result_list.append({
        'Stat1': sum(sub_data['Var1']),
        'Stat2': min(sub_data['Var1']),
    })
result = pd.DataFrame(result_list)
```

## 3. 缺失处理相关

### 3.1. 读入时缺失指定

通过`na_values`参数可以将多种数值替换为缺失:

```python
df = pd.read_table(FILE, sep="\t", header=0, na_values=['.', ''])
```

### 3.2. 输出时缺失指定

输出时可通过参数指定缺失值用什么来填充：

```python
df.to_csv(FILE, sep="\t", header=0, na_values=['.', ''])
```

### 3.3 是否缺失判定

pandas内如果产生缺失, 会使用numpy内置的`nan`对象, 这个对象在判定上和R中存在的`NA`或`NULL`以及Python内置的`None`都不太一样, 因为从类型上`nan`居然是`float`... 通过`A==nan`这种逻辑判断是得不到`True`的结果的... 原因我暂时没有找到. 不过要判断是否缺失的话, 最好还是用pandas下面的`isnull`方法. 其他模块像`math`, `numpy`虽然也有内置缺失判断函数, 但是这些函数有个致命的缺陷...他们接受的参数类型只能是`float`... 所以如果用他们的我还得自己加一层类型判断, 着实有点麻烦, 所以还是用pandas下的`isnull`最方便.

## 4. 最常用的输出参数

因为个人习惯性输出tsv文件, 因此用pandas输出文件的话最经常这么写:

```python
df.to_csv(FILE_PATH, sep="\t", index=0, na_rep=".")
```

需要注意的是, 如果输出的结果表是数据透视表, 很大概率是要将`index`输出的, 这种时候`index=0`就要去掉了

## 5. 其他可能需要的操作

### 5.1 列类型转换

```python
df.astype({'col1': 'int32'})
```
