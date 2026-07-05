---
title: 把 Claude Science 装进 Docker 容器
date: 2026-07-05 23:00:00
tags:
  - Docker
  - Claude Science
  - API
  - AI
categories: AI
---

前几天折腾了一个小项目：把 Claude Science 的 Linux 版本装进 Docker 容器，再配上 API Bridge，让它能通过第三方 API 后端运行，不需要直接持有 Anthropic 的 API Key。

过程比预想的曲折，但也挺有意思。记录一下。

<!-- more -->

## 项目架构

整个架构其实不复杂：

```
┌─────────────────────────────────────────┐
│              Docker Container            │
│                                          │
│  ┌──────────────┐    ┌────────────────┐  │
│  │  API Bridge  │◄───│ Claude Science  │  │
│  │  (:9876)     │    │  (:9981)        │  │
│  └──────┬───────┘    └────────────────┘  │
│         │                                │
└─────────┼────────────────────────────────┘
          │
     ┌────┴────────────┐
     │ Third-Party API │
     │ (SiliconFlow,   │
     │  Moonshot, etc) │
     └─────────────────┘
```

Claude Science 通过 operon 守护进程运行，它会向 [API Bridge](https://github.com/SilenWang/Claude-Science-Container) 发送 Anthropic 格式的请求。API Bridge 负责将请求转换为 OpenAI 兼容格式，再转发给后端的第三方 API。

## 开发过程中的坑

### 1. Claude Science 的访问地址限制

Claude Science中给了 `--host` 参数可以指定能访问的地址，但是实际上，并不能通过 localhost 以外的ip进入软件界面，所以要么端口转发作为本地使用，要么就得穿透服务和这个容器一台机器，否则根本无法使用...

### 2. https链接限制

原本我是按照一个服务一个容器的原则，将api bridge和claude science分开放，设置好了发现 Claude Science 发请求的url限制https或本地，只好又回退到一个容器...

### 3. protocol 转换的限制

API Bridge 在做 Anthropic -> OpenAI 格式转换时，会丢失 Anthropic 特有的一些功能。最明显的就是 **tool calls** 和 **web_search**。OpenAI 兼容的接口没有对应的实现，所以这些功能就用不了。

不过 DeepSeek 提供了一个 Anthropic 兼容的 direct API 端点，如果走这个通道就可以绕过转换，所有功能都完整保留。但代价是模型选择范围就窄了——只能选 DeepSeek 官方的模型API。

### 4. MCP 目录连接器的问题

进入 Claude Science 的设置界面时，会提示 "Directory connectors unavailable"。一开始以为是配置有问题，查了一圈日志发现是 `claudeAiFetch: 401 and refresh failed`——它试图去 claude.ai 同步 MCP 服务器列表，但因为没有登录 Anthropic 账号，所以失败了。

不过实际上不影响核心功能。Claude Science 内置的科学 MCP 服务器（biorxiv、chembl 这些）依然能正常连接和使用。类似于"报了个错但其实不影响使用"的情况。

这让我想起以前用 snakemake 时，有些规则的 shadow 功能在特定版本会报 warning 但不影响执行。当时觉得这种半残废的功能很烦人，但也学会了怎么快速判断哪些 warning 可以忽略。这次遇到 MCP 目录的报错，扫了一眼日志就知道是 auth 的问题而非功能性问题，也就不纠结了。

### 5. 任务莫名终止

用 DeepSeek Flash 跑任务时，有时任务会中途莫名其妙终止，但 DeepSeek Pro 就好很多。这个问题定位起来比较困难，因为既不是配置问题也不是网络问题，更像是模型自身的行为差异。只能认为是替换的模型智力还不到了...

## 项目地址

代码放在 GitHub 上：[Claude-Science-Container](https://github.com/SilenWang/Claude-Science-Container)

欢迎有兴趣的朋友试试，有问题可以提 issue。

## 小吐槽

之前写过一篇博客吐槽 API 互不兼容的问题，说 OpenAI 和 Anthropic 的 API 接口不兼容，切换模型要改一堆代码。当时还感叹，这跟打印机厂商搞耗材垄断一个套路。

结果没想到，这个吐槽很快就变成了实际要处理的问题——Claude Science 作为 Anthropic 的产品，它的 API 调用自然是 Anthropic 风格的。但我想用的是第三方的 API 后端（OpenCode 这些），它们只支持 OpenAI 兼容格式。也就是说，中间需要一个翻译层，把 Anthropic 风格的请求转成 OpenAI 兼容的格式。还好社区大神们给力，一夜就有了方案，我只用做个调包侠，就能解决问题了...