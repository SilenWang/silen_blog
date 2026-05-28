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

最近帮朋友做了个小工具，需要在浏览器里运行 Python 代码进行数据处理，而且用户环境可能没有网络。这个需求放在几年前几乎不可能实现——要么你得搭个后端服务器跑 Python，要么就得用纯 JS 重写整套逻辑。但有了 Pyodide + PWA 这套组合拳，事情就变得可行了。

<!-- more -->

## 先说说背景

朋友的需求是这样的：有一个比较特定的数据处理流程，代码是用 Python 写的，依赖了一些常用的科学计算库（pandas、numpy 这类）。用户需要在浏览器里上传文件、处理、下载结果，而且很多时候用户是在没有网络的内网环境或者飞行模式下使用。

这个场景其实挺典型的。很多生物信息、数据分析领域的小工具，核心逻辑都是 Python 写的，但分发方式如果不是让大家装 Python 环境跑脚本，就是搭个 Web 服务。前者对非技术用户不友好，后者需要服务器成本和运维。

## PWA 构建曾经的门槛

PWA（Progressive Web App）的概念已经出来好几年了，核心能力就是通过 Service Worker 实现离线访问、后台同步等功能。但过去想做一个实用的 PWA 工具，面临两个比较大的问题：

### 1. 业务逻辑必须用 JS 重写

如果你已有的核心逻辑是 Python 写的，想在浏览器里跑，基本上只有两条路：

- **重写为 JavaScript**：听起来简单，但 Python 的科学计算生态（pandas、numpy、scipy 等）在 JS 里没有完全对等的替代品。即使有类似功能的库，行为细节的差异也够你调试好一阵子。
- **架设后端服务**：前端只负责展示，后端跑 Python 逻辑，通过 API 通信。这个方法可行，但引入了服务器成本、网络依赖和运维复杂度。而且如果用户需要离线使用，这招就行不通了。

### 2. HTTPS 的硬性要求

PWA 的核心特性（Service Worker）要求在 HTTPS 环境下才能注册（localhost 除外）。这就意味着，如果你想做一个 PWA 工具给人用，你得先解决 HTTPS 的问题。

对于有固定域名的生产环境，申请 SSL 证书不是问题。但如果你只是想做个临时工具分享给朋友或同事用，就有点尴尬了：

- 为了一个临时工具去申请域名和证书，太折腾
- 用自签名证书，浏览器会报不安全
- 用 ngrok 这类工具可以快速暴露 HTTPS 端口，但依赖第三方服务

实际上，在开发测试阶段，最好的方案就是直接用 `localhost`（浏览器对 localhost 豁免了 HTTPS 要求），然后用 `--host` 参数让同局域网的人访问。但这种方式只适用于局域网，而且仅限于开发测试。

## Pyodide：在浏览器里跑 Python

[Pyodide](https://pyodide.org/) 是一个将 CPython 编译到 WebAssembly 的项目，让你可以在浏览器里直接运行 Python 代码。它不仅仅是把 Python 解释器搬到了浏览器里，还通过 Emscripten 做了大量的适配工作，让 numpy、pandas、scipy、matplotlib 这些常用的科学计算库也能在浏览器里跑。

使用方式很简单，在 HTML 里引入：

```html
<script src="https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js"></script>
<script>
  async function main() {
    let pyodide = await loadPyodide();
    await pyodide.loadPackage(['numpy', 'pandas']);
    // 现在可以直接在浏览器里跑 Python 了
    let result = await pyodide.runPython(`
      import numpy as np
      import pandas as pd
      data = pd.DataFrame({'a': [1, 2, 3], 'b': [4, 5, 6]})
      data.describe()
    `);
    console.log(result);
  }
  main();
</script>
```

当然，实际使用中你不会直接写大段字符串来传 Python 代码。更合理的做法是把 Python 脚本作为独立文件，通过 Pyodide 提供的文件系统 API 加载执行。

```javascript
// 将 Python 文件写入 Pyodide 的虚拟文件系统
pyodide.FS.writeFile('/main.py', pythonCode);

// 执行
let result = pyodide.runPython(`
  import sys
  sys.path.insert(0, '/')
  import main
  main.run()
`);
```

## 结合 PWA 实现离线能力

有了 Pyodide 解决了"浏览器里跑 Python"的问题，再结合 PWA 的离线能力，就能做到完全不依赖网络的完整工具链了。

整体架构大概是这样：

1. **核心逻辑**：Python 脚本，通过 Pyodide 在浏览器端执行
2. **前端界面**：HTML + CSS + JS，提供文件上传、参数配置、结果展示
3. **Service Worker**：缓存所有静态资源（HTML、JS、CSS、Pyodide 的 wasm 文件、Python 包等），确保离线可用
4. **IndexedDB**：如果需要持久化用户数据或者缓存处理结果，可以存到 IndexedDB 里

这里有个关键点：Pyodide 的 wasm 文件和 Python 包（.whl 或 .data 文件）体积不小，首次加载时需要网络。但一旦通过 Service Worker 缓存下来，后续即使完全离线也能正常使用。

```javascript
// Service Worker 中的缓存策略
const CACHE_NAME = 'pyodide-app-v1';
const PRECACHE_URLS = [
  '/',
  '/index.html',
  '/app.js',
  '/style.css',
  '/main.py',
  // Pyodide 相关资源
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.js',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.wasm',
  'https://cdn.jsdelivr.net/pyodide/v0.25.0/full/pyodide.asm.data',
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_NAME).then(cache => cache.addAll(PRECACHE_URLS))
  );
});
```

## 临时 HTTPS 访问的解决方案

前面提到 PWA 需要 HTTPS，而临时工具有时候不方便搞域名证书。这里有几个我试过的方案：

### 方案一：localhost + 局域网共享

如果是开发测试或者小范围使用，直接在本地起一个 HTTPS Server（用自签名证书或者 mkcert 生成可信的本地证书），然后用 `--host 0.0.0.0` 监听所有网络接口，局域网内的设备就可以通过 `https://192.168.x.x:port` 访问了。

```bash
# 使用 mkcert 生成本地可信证书
mkcert -install
mkcert localhost 192.168.1.100 127.0.0.1

# 用 node 的 https 模块或者 caddy 起服务
caddy file-server --domain localhost
```

mkcert 生成的证书会被系统信任，浏览器不会报不安全，体验上接近生产环境。

### 方案二：Cloudflare Tunnel

如果你需要一个公网可访问的 HTTPS 地址做演示（比如给不同城市的朋友测试），Cloudflare Tunnel 是个好选择。它会建立一个到 Cloudflare 边缘节点的隧道，自动提供 HTTPS。

```bash
cloudflared tunnel --url http://localhost:8080
```

执行后会得到一个 `https://xxx.trycloudflare.com` 的地址，直接在浏览器打开就行，HTTPS 和 PWA 都能正常工作。不需要域名，不需要配置证书，用完关掉就行。

### 方案三：GitHub Pages

如果工具本身是纯前端的（不涉及后端 API），那也可以直接部署到 GitHub Pages。它自动提供 HTTPS，而且有固定的 URL，适合长期可用的工具。不过缺点是更新需要走 git 流程，不如本地开发灵活。

## 实际效果和体验

经过这样一套组合，最终的效果是：用户第一次访问时需要联网加载资源（主要是 Pyodide 的 wasm 文件和 Python 包），之后即使完全断网，刷新页面也能正常使用所有功能。

处理速度方面，Pyodide 通过 WebAssembly 运行 Python，性能虽然比不上原生 CPython，但对于大多数数据处理任务来说已经够用。尤其是在数据量不大的场景下（几 MB 到几十 MB 的表格数据），处理时间基本在秒级，用户感知不到明显延迟。

## 总结

Pyodide + PWA 这套组合，让"在浏览器里离线运行 Python 代码"这件事从一个听起来很酷但不太现实的想法，变成了一个真正可落地的方案。

它的核心优势在于：

- **无需后端**：纯前端架构，省去了服务器成本和运维
- **离线可用**：Service Worker 缓存后，完全断网也能工作
- **复用 Python 生态**：不需要把已有的 Python 代码重写成 JS，numpy、pandas 这些库直接用
- **零部署成本**：静态文件往哪一扔都能跑，GitHub Pages、本地服务器、甚至 U 盘里打开都行

当然它也有局限：不适合计算量大的任务（Pyodide 的性能天花板在那里），不适合需要调用系统 API 的场景（浏览器沙箱限制），首次加载的资源体积也不小。

但对于"有一个 Python 写的工具，想让非技术用户在浏览器里离线使用"这个场景来说，目前可能是最好的方案了。
