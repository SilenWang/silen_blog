---
title: Python命令行解析
categories: Script
date: 2018-10-22 12:12:50
tags: ['python', 'argparse', '参数解析']
---

一直都在用`argparse`写命令行解析, 但是直到前段时间才发现原来它可以用来写子命令
<!-- more -->

## 一般参数解析

一般的参数解析很简单了, 示例化一个`argparse`, 然后调用对象的添加参数方法就可以

```python
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(
    prog='program',
    description='Some description',
    formatter_class=RawTextHelpFormatter
    )
parser.add_argument('-s',
                    help="description for args",
                    required=False,
                    type=str)
```

## 子命令创建

子命令需要在主参数解析对象(上面是`parser`)的基础上首先创建子解析对象,然后调用这个对象的方法添加子命令

```python
sub_parser = parser.add_subparsers(title='sub-commands',
                                    description='Some description',
                                    help='Some help description')
sub_parser_a = sub_parser.add_parser('a', help='Some help description')
sub_parser_a.add_argument('-s', 
                        help="description for args",
                        type=str)
```

这样就可以实现在主程序下定义子程序了, `argparse`还提供了一个方法, 在特定子命令被调用时, 指定运行某函数, 这个函数会接收到所有子命令下参数值组成的一个命名空间, 可以进行参数解析后开始执行子程序

```python
sub_parser_a.set_defaults(func=subFun)

def countFun(args):
    args = vars(args)
    print(args)
```

## 其他

`argparse`还可以实现多重子程序以及子程序参数继承主程序参数, 这个还未实际使用过, 后面再更新~