---
title: Let Custom Objects Use the 'with' Statement
categories: Coding
date: 2019-11-14 18:56:24
tags: ['Python']
---

The `with open` statement is a particularly convenient feature in Python, allowing for automatic closure of objects when operations are completed. For custom objects, it's also possible to support the `with` syntax.
<!-- 摘要部分 -->
<!-- more -->

To achieve this, define the `__enter__` and `__exit__` methods in your class:

```python
class MyClass(object):
    def __init__(self):
        pass
    def __enter__(self):
        print("Start")
    def __exit__(self, exc_type, exc_val, exc_tb):
        print("End")
```

With this setup, when using the `with` statement, it will first execute the content within `__enter__`, and automatically execute the content within `__exit__` at the end of the code block.

```
>>> with MyClass() as tClass:
...     print("test")
...
Start
test
End
```
```
