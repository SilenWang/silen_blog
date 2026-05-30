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

最近帮朋友做了个小工具：处理输入文件中的序列，生成指定格式的 XML。朋友不习惯命令行，需要图形界面；且输入数据有保密要求，计算逻辑必须在本地完成。按常规该开发桌面软件，但考虑到后续维护更新，我想试试可离线的 PWA Webapp。查了下，Pyodide 中已有 pandas 的 wasm，运行库无需额外工作——那就这么干！
 
<!-- more -->

## PWA 构建曾经的门槛

PWA（Progressive Web App）概念提出已久，核心能力是通过 Service Worker 实现离线访问、后台同步等功能。但过去想做一个离线 PWA 工具，对我来说相当有难度——我不会从头写 JS 或 TS。

## Pyodide：在浏览器里跑 Python

[Pyodide](https://pyodide.org/) 是一个将 CPython 编译到 WebAssembly 的项目，让我们可以在浏览器里直接运行 Python 代码。它不仅把 Python 解释器搬到了浏览器中，还通过 Emscripten 做了大量适配，让 numpy、pandas、scipy、matplotlib 这些科学计算库也能在浏览器里跑——它们底层并非 Python，光有解释器无法运行。

## 结合 PWA 实现离线能力

Pyodide 解决了"浏览器里跑 Python"的问题，再结合 PWA 的离线能力，就能做到完全不依赖网络的完整工具链。

整体架构大概是这样：

1. **核心逻辑**：Python 脚本，通过 Pyodide 在浏览器端执行
2. **前端界面**：HTML + CSS + JS，提供文件上传、参数配置、结果展示
3. **Service Worker**：缓存所有静态资源（HTML、JS、CSS、Pyodide 的 wasm 文件、Python 包等），确保离线可用

Pyodide 的 wasm 文件和 Python 包（.whl 或 .data）体积不小，首次加载需要网络。但一旦通过 Service Worker 缓存下来，后续即使完全离线也能正常使用。

## AI 帮我完成主要编码工作

有了前面的框架，只能说想法行得通。但换作过去，写这个程序依然费劲——需要多语言编码，且用 Python 框架做 Web 界面不划算（会缓存更多 Python 库，还可能引入难以编译成 wasm 的库）。

但现在有 DeepSeek 了。这个程序需要的界面不过两个按钮加一个日志框，完全交给它帮我完成就行。

## 实际效果和体验

最终效果：用户首次访问时需联网加载资源（主要是 Pyodide 的 wasm 文件和 Python 包），之后即使完全断网，刷新页面也能正常使用所有功能。

处理速度方面，通过 WebAssembly 运行 Python 性能不及原生 CPython，但这个工具处理的数据不到 1 万条，无需多高的性能。

以后我也可以用这套组合开发更多小工具。