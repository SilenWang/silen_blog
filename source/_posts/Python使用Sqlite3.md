---
title: Python使用Sqlite3
categories: Script
date: 2018-08-26 23:54:24
tags: ['python', 'sqlite']
---

为了以后操作使用Sqlite3做准备, 上周上班的时候特意试了一下用操作数据库的数据, 为了方便直接用了sqlite3, 没想到sqlite3并不是所有语句都支持, 并且sqlite自由的点命令是无法被python的sqlite模块调用的....

<!-- more -->

## 简单的使用方式:

```python
import sqlite3   # 导入模块
conn = sqlite3.connect('test.db')   # 连接数据库(如果不存在的话会被创建)
cursor = conn.cursor()    # 从数据库对象创建指针
cursor.execute('sql cmd')   # 向数据库传入sql语句, 注意这里的操作不是立刻保存/同步到数据库中
cursor.close()   # 关闭指针
conn.commit()   # 确认更改
conn.close()   # 关闭连接
```

- 使用中的注意事项
    * 只支持标准的SQL语句, 由于使用的是sqlite, 其不支持的SQL语法一概无法使用, 也不支持使用sqlite特有的点命令, 因为这些命令本质上是在shell内执行的, 不是SQL语句
    * 如果一定要使用点命令, 则需要使用os模块或其他可以调用shell的模块, 本质上是在执行shell命令
    * 不是所有的操作都一定要指针对象完成, 有部分可以使用数据库对象, 这样可以少写两句, 具体见[菜鸟教程网站的介绍](http://www.runoob.com/sqlite/sqlite-python.html)
    * 注意一定要确认更改(`obj.commit() `), 否则数据库内的内容是不会变化的
    * 尝试使用了with的方式进行连接, 但是在with语句的范围外是不会自动断开连接的(`obj.close()`), 所以还是得老实写

## 调用点命令的方式

```python
import os, textwrap

# 这里是使用textwrap进行了shell命令的缩进调整, 并不必要, 只是输出的时候方便看
def cmd_gen(cmd_str):
    cmd_str = textwrap.dedent(cmd_str)
    cmd_str = cmd_str.strip()
    return cmd_str


db_file = 'test.db'
index_file = 'index.txt'

# 生成命令, 这个写法是我从网上搜回来的, 是找到的唯一有效的写法, 具体意义日后再看...
cmd = cmd_gen('''
    sqlite3 {db_file} << EOF
    .separator "\\t"
    .import {index_file} idx_tab
    EOF
    '''.format(db_file=db_file, index_file=index_file))
os.system(cmd)
```