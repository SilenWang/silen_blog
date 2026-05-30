---
title: 利用 Pyodide 构建可离线使用的 PWA 应用
date: 2026-05-28 07:30:00
tags:
  - Pyodide
  - PWA
  - WebAssembly
  - Python
  - Service Worker
  - JavaScript
categories: Coding
---

最近帮朋友做了个小工具，主要功能是处理输入文件中的序列，生成指定格式的XML文件。朋友并不习惯使用命令行，因此需要图形界面，然后被处理的输入有保密要求，因此计算逻辑需要在本地设备完成。在这个要求下，一般就应该开发一个桌面软件了，但是考虑到后续维护更新方便，我想试试写可以离线的PWA Webapp，查了一下，Pyodide中已经有pandas的wasm了，因此运行库上我不需要做什么了，那就这么干！
 
<!-- more -->

## PWA 构建曾经的门槛

PWA（Progressive Web App）的概念已经出来非常长时间了，核心能力就是通过 Service Worker 实现离线访问、后台同步等功能。但过去想做一个能离线的 PWA 工具，对我来说相当有难度，毕竟，我并不会从头写Js或Ts。

## Pyodide：在浏览器里跑 Python

[Pyodide](https://pyodide.org/) 是一个将 CPython 编译到 WebAssembly 的项目，让我们可以在浏览器里直接运行 Python 代码。它不仅仅是把 Python 解释器搬到了浏览器中，还通过 Emscripten 做了大量的适配工作，让 numpy、pandas、scipy、matplotlib 这些常用的科学计算库也能在浏览器里跑（他们底层并不是Python，光有Python解释器无法运行）。

## 结合 PWA 实现离线能力

有了 Pyodide 解决了"浏览器里跑 Python"的问题，再结合 PWA 的离线能力，就能做到完全不依赖网络的完整工具链了。

整体架构大概是这样：

1. **核心逻辑**：Python 脚本，通过 Pyodide 在浏览器端执行
2. **前端界面**：HTML + CSS + JS，提供文件上传、参数配置、结果展示
3. **Service Worker**：缓存所有静态资源（HTML、JS、CSS、Pyodide 的 wasm 文件、Python 包等），确保离线可用

这里有个关键点：Pyodide 的 wasm 文件和 Python 包（.whl 或 .data 文件）体积不小，首次加载时需要网络。但一旦通过 Service Worker 缓存下来，后续即使完全离线也能正常使用。

## AI 帮我完成主要编码工作

有了前面的框架，只能说我的想法行得通，但是如果是过去，写这么一个程序依然是很费劲的，毕竟依然需要多语言编码（再用Python下的框架去实现Web界面是不划算的，会需要缓存更多的Python库，同时也可能引入难以编译成wasm的库）。

但是现在有 DeepSeek 了，介于这个程序需要的界面，两个按钮一个日志框就完事了，完全可以交给Deepseek帮我完成。

## 实际效果和体验

经过这样一套组合，最终的效果是：用户第一次访问时需要联网加载资源（主要是 Pyodide 的 wasm 文件和 Python 包），之后即使完全断网，刷新页面也能正常使用所有功能。

处理速度方面，虽然通过 WebAssembly 运行 Python，性能远比不上原生 CPython，但这个工具需要处理的数据低于1万条，并不需要多高的效能。

我想之后我也可以用这套组合，开发更多的小工具。