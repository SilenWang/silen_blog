---
title: 还是搞不懂_为什么taipy的on_action不嵌套一下就报错
categories: Coding
date: 2025-12-17 21:45:09
tags: ['taipy', 'callback']
---

我使用过的所有应用快速开发框架，清一色的都是将元素动作绑定到某个Python函数，从而触发信息的更新或者改变，所以应该还算比较有经验。但是Taipy这个的属实让我有点摸不着头脑，时不时的就会报错...

<!-- more -->


在 Taipy 中，当我们为 `on_action`（或 `on_change`）属性直接传递一个函数对象时，会很随机的遇到类似“函数不合法（function not valid）”的错误。例如，以下写法可能会报错：

```python
tgb.table(
    ...
    on_action=update_prod_link,  # 直接传递函数名
    ...
)
```

但是如果用匿名函数，就不会出现前述的问题，一切正常：

```python
tgb.table(
    ...
    on_action=lambda s, v, p: update_prod_link(s, v, p),
    ...
)
```

由于官方的例子都比较简单，都是写的匿名函数，Issue区也没有看到类似的问题反馈，目前阶段我也就只能这么脱裤子放屁了...