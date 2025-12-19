---
title: 我真是用什么都会踩坑--Devcontiner使用的若干问题
categories: Others
date: 2025-12-19 22:27:38
tags: ['devcontainer', 'devpod', 'docker', 'git']
---

我只是想设置一个devcontainer环境，更高效的完成维护公司官网的工作，没想到一件事上能踩三个坑... 

<!-- more -->

## 背景

[Devcontainer](https://code.visualstudio.com/docs/devcontainers/containers) 是 Visual Studio Code 推出的一项功能，它允许我们利用 Docker 容器快速搭建与项目匹配的开发环境。这种“代码即环境”的方式可以极大地统一团队协作体验，减少“在我机器上能跑”的尴尬。我之前写的devpod就是基于这种格式的开源解决方案。

然而，理想很丰满，现实却总会在细节上给你使绊子。在最近一次为公司官网项目配置 Devcontainer 时，我接连遇到了三个意想不到的问题，每个都耗费了不少时间才找到解决方案。下面就把这些问题和解决方案记录下来，供大家参考。

## 坑一：Git 子模块（subtree）的缺失

### 问题描述
微软官方提供的Devcontainer基础镜像有 Debian 和 Ubuntu，两者都自带的 Git，但是，版本都比较基础。不含 subtree 子命令... 嗯，我也是第一次知道，git 的功能太多，以至于一些比较新的命令，做成了可选项... 

但这都不是问题，Ubuntu 和 Debian 源中的 Git 应该都是功能完整的，因为我在 Devcontainer 之外从来没有收到过 subtree 命令不存在的提示。然后问题就出现了，即使通过 `apt-get install git-man` 或者 `git-all`，也无法直接获得 `subtree` 功能。为什么呢？其实很简单，Devcontainer 内自带的git，并不来自软件源，它位于 `/usr/local/bin`，而系统安装的 git 会在 `/usr/bin`。而环境变量 `Path` 中，`/usr/local/bin` 优先级更高， 所以自行安装的git，实际上不会被使用到，于是会一致找不到 `subtree` ...

### 解决方案
根据前面的问题，其实通过 `apt` 安装系统带的 git，然后把镜像里的干掉就是：

```dockerfile
FROM mcr.microsoft.com/devcontainers/base:ubuntu

RUN apt-get update && \
   add-apt-repository -y ppa:git-core/ppa && \
   apt-get update && \
   apt-get install -y git vim && \
   # 记得干掉容器里原来的git
   rm -r /usr/local/bin/git* && \
   apt-get autoremove -y && \
   apt-get clean -y
```

## 坑二：Windows 平台下的换行符自动转换

### 问题描述
当你在 Windows 主机上使用 Docker Desktop 运行 Devcontainer 时，Git 会根据 `core.autocrlf` 的默认设置（通常为 `true`）自动将行结束符 LF 转换为 CRLF。这本身是 Windows 开发环境中的常见行为，但在容器内部，文件系统是 Linux 风格的，因此 Git 会在检出文件时执行转换，导致几乎所有文件都被标记为“已修改”。

具体现象是：在容器内执行 `git status`，会发现大量文件显示为修改（即使你什么都没做）。这种干扰不仅让人困惑，还可能影响后续的提交操作。

### 解决方案
Emmmm，我其实并没有解决这个问题，因为进行 git 的配置并不能阻止 vscode 执行这个转换... 这个问题只发生在 Windows 平台下使用 Docker Desktop 的情况，反正用 Docker Desktop 编译个网站都能导致桌面卡死，直接不用就是了... 老老实实用服务器作为 SSH Provider...

## 坑三：Devcontainer Desktop 与远程 SSH Provider 的掉线问题

### 问题描述
在尝试通过 Devpod Desktop 连接到远程 SSH Provider时，会发现连接极其不稳定，经常在几分钟内就会自动断开，且重连过程十分缓慢，连上不一会又掉。不足。由于 Devpod Desktop 也不给什么运行日志，我也不知道是啥毛病...

### 解决方案
简单粗暴，Devpod Desktop 只用来进行 Provider 和工作区的配置，发起连接由 DevPod CLI 完成就好...

## 小结
好消息，问题最终都解决了；坏消息，本周又少睡了七八个钟，愿我不要暴毙在2025...