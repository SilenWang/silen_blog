---
title: 搭一个用rendercv生成简历的容器
categories: Coding
date: 2026-01-06 21:17:55
tags:
  - RenderCV
  - 容器
  - 开发环境
  - Pixi
  - Aider
---

准备简历说起来是不是什么难事，但实际做起来，还是会比较耗时的。现在不少岗位会对细分领域有一些只是或者项目的要求，求职者这么多，一份没有针对岗位要求突出重点的简历，也许就被用人单位忽略掉了。所以想要达到更好的效果，最好能对每个岗位针对性的优化... 这一点AI应该是擅长的，但是至少我现在没找到比较好的免费工具。

刚好最近看到了[RenderCV](https://github.com/rendercv/rendercv)，这是一个通过YAML配置文件生成简历的工具，它可以实现配置文件即简历，将它与[Aider](https://aider.chat/)这种编程助手结合，基本就能构成一个快速准备简历的工作环境了。

<!-- more -->

## 容器的构成组件

容器的组件很简单：

- Pixi：用于安装依赖，即RenderCV和Aider
- Aider：用于生成大致的简历配置文件，同时也可以直接用来翻译
- RenderCV：从配置文件生成简历
- PDF预览插件：直接在VSCode中就能查看生成的PDF简历，配合RenderCV的`--watch`参数，也可以实现修改的实时预览

## 启动容器

- **使用Github Codespace**：最方便的方式当然是直接白嫖Gihub的Codespace，fork [我的容器项目后](https://github.com/SilenWang/RenderCV_Pod)，就可以在Web版的VSCode里准备简历了。简历生成后也可以直接从容器里下载到本地。

- **使用Devpod**：使用Devpod也简单，有任意Provider后，运行 `devpod up https://github.com/SilenWang/RenderCV_Pod`就能开始写了

## Aider生成简历配置的注意事项

目前我配合DeepSeek的API进行简历生成，可能由于RenderCV项目本身比较新，DeepSeek显然是不知道RenderCV的配置文件应该是什么样的，它能给出看上去可以的文件，但是其中Section部分该怎么写完全是错误的，还是需要自行看[RenderCV的文档](https://docs.rendercv.com/)来修正，当然RenderCV的设计就是为了简单快速，所以这些调整10分钟左右就能学会，到不是太大问题。

另外，准备一个示例完整的例子给AI以后，基本就能避免前述问题，之后有空我再在项目中补充。