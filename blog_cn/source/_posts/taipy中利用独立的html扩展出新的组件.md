---
title: taipy中利用独立的html扩展出新的组件
categories: Coding
date: 2025-12-17 21:42:31
tags: ['taipy', 'html']
---

我之前就在疑惑，Taipy看上去是在活跃维护的，用户看上去应该也不少（一堆用户提的Issue），但是它的组件库并不能算多么丰富，甚至有些基本的功能都小bug一堆（比如我之前发现的明暗切换按钮bug）。然后我无意间看到了官网上嵌入第三方内容的教程，我突然理解了，虽然它组组件少，但是他可以缝啊！

<!-- more -->

## 为什么需要嵌入第三方组件

Taipy 作为一个以 Python 为中心的 Web 应用框架，其内置的视觉元素（Visual Elements）已经能够覆盖很多常见的交互场景。但在实际项目中，总能轻松的从甲方那收到刚好在现有组件外的要求，这时候，我们只能自己想办法将这些内容“缝”进 Taipy 的页面里。

Taipy固然提供扩展当前Python组件的说明，但是为了一个用过即丢的小项目，实在没必要先写html/js，然后准备好接口绑定到Taipy，然后再做测试、编译、重新安装... 如果真这么干，感觉完成组件前我就被老板开了... 

幸运的是，Taipy 还提供了另外一种“缝”的机制，允许我们将任意 HTML + css + js 编写的页面直接在页面中渲染（其实就是弄了个iframe 套进去）。这个机制的核心就是 `Gui.register_content_provider()` 函数。

## 核心机制：register_content_provider

`Gui.register_content_provider(type, provider_func)` 接受两个参数：

- `type`：需要是一个Python类型（例如 `folium.folium.Map`）
- `provider_func`：一个回调函数，它接收一个该类型的对象，并返回一个 **bytes** 对象，这个 bytes 对象应当是该对象对应的 HTML 内容。

当 Taipy 在页面中遇到一个属于该类型的变量，并且这个变量被放在 `part` 组件的 `content` 属性中时，Taipy 就会自动调用你注册的 `provider_func`，将返回的 HTML 以 iframe 的方式嵌入页面中。

这样一来，只要能最终得到 HTML 文件 的内容 Taipy 就能进行渲染，显示在页面上。

官方给的例子是[Folium地图对象](https://docs.taipy.io/en/latest/tutorials/articles/3rd_party_components/)，这个库与plotly类似，也是将多种地图的js库进行封装后，提供python api，完成绘制后，可以直接输出html页面，因此可以将生成的html内容直接嵌入到taipy中。

## AI 带来了更强的扩展性

如果是在两年前，Taipy的方案其实多少有些鸡肋，因为扩展功能最终还是要老老实实的编写Web那一套，违背通过Python快速开发的初衷。但是 LLM 当道，情况发生了极大的改变，现在各种 LLM 可以在极短的时间内生成一个可用的 Web 小组件，然后只需要稍微学习 [Jinja2模板](https://docs.jinkan.org/docs/jinja2/) ,就能快速的将数据和Web小组件连接起来，扩展Taipy的功能。虽然这不适用于要求功能特别复杂的情况，但是开发个原型程序绝对够用了。
