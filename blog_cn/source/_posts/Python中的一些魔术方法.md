---
title: Python中的一些魔术方法
categories: Script
date: 2019-11-25 23:34:06
tags: ['Python', '魔术方法', '面向对象', '特殊方法']
---

原来之前使用的`__enter__`和`__exit__`就是以前耳闻的魔术方法, 这次我的程序开发因为以面向对象为主, 为了实现一些python中好用的功能所以涉及到了更多方法, 这里做个小笔记
<!-- 摘要部分 -->
<!-- more -->

### `__len__`

定义之后, 可以对实例使用`len()`函数.

```python
def __len__(self):
    '''
    定义长度
    '''
    return len(self.Seq)
```

### `__getitem__`

定以后, 可以对实例进行`[]`索引和切片操作

```python
def __getitem__(self, key):
    '''
    添加切片索引
    '''
    return self.d[key]
```

### `__str__`

定义后, 使用`print()`打印实例时会打印这个方法内返回的内容

```python
def __str__(self):
    return str(self.__dict__)
```

### `__iter__`

定以后, 可对对象进行迭代操作, 比如`for i in object`这种. 需要注意的是, 这个方法必须要返回迭代器(可迭代对象不行).

```python
def __iter__(self):
    return (mut for mut in self.muts)
```


### `__eq__`和`__hash__`

这两个放一起是因为把实例放在`set()`中, 想要自动去重时, 这两个都要有, 否则无法达到目的效果. 其中`__eq__`是使对象可被比较是否相同(`==`和`!=`, `is`好像不可以). 注意这个方法里的`other`好像是固定的, 不要改动.

```python
def __eq__(self, other):
    '''
    保证利用set特性去重
    '''
    return self.muts == other.muts
```

`__hash__`则本质上是对象可被哈希化, 哈希化时的哈希值会是这个方法给的哈希值

```python
def __hash__(self):
    '''
    保证利用set特性去重
    '''
    return hash(self.muts)
```

暂时就这么多, 以后有用到新的继续水~
