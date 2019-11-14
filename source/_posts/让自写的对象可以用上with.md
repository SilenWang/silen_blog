---
title: 让自写的对象可以用上with
categories: Script
date: 2019-11-14 18:56:24
tags: ['Python']
---

`with open`是python中特别方便的一种写法, 对于支持的函数来说, 有了`with`这个语句可以很方便的在操作完成时自动关闭对象. 对于自己编写的对象其实也可以实现`with`语法的支持.
<!-- 摘要部分 -->
<!-- more -->

具体的做法就是在自己的`class`中定义`__enter__`和`__exit__`两个方法:

```python
class MyClass(object):
    def __init__(self):
        pass
    def __enter__(self):
        print("Start")
    def __exit__(self, exc_type, exc_val, exc_tb):
        print("End")
```

如此一来, 在使用`with`时就会先执行`__enter__`内的内容, 在代码块结束的时候则自动执行`__exit__`的内容了.

```
>>> with MyClass() as tClass:
...     print("test")
...
Start
test
End
```