---
title: 我的 channel 的二三事
date: 2026-05-08 07:00:00
tags:
  - channel
  - pixi
  - GitHub Actions
categories: Others
---

我在 [prefix.dev](https://prefix.dev) 上维护了一个 软件包 channel，用来分发一些自己打包的软件。由于最近都在忙工作相关的事情，有一段时间没有管它了，今天一看，Action 居然莫名其妙不运行了... 于是了解了下原因，顺便记录几个其他的问题。

<!-- more -->

## 配额限制

prefix.dev 上的 open channel 虽然给的非常慷慨，但并非资源无限。每个 channel 的存储空间和包数量都有配额限制。目前我的channel只有5个包，却也消耗了10G的空间，主要还是opencode更新的非常频繁，一个小版本下都能有几十个子版本，虽然一次打包也就几十m，但架不住版本真的多。

目前我得手动检查 channel 里的旧包，手动删除那些不再需要或者已经被新版本替代的包。这个过程完全是手工操作，比较繁琐，但暂时没有找到更好的自动化方案。

## GitHub Actions 的定时任务会被自动关闭

我依赖 GitHub Actions 来自动构建和更新 channel 中的包。但 GitHub 有一个策略：如果一个仓库的 Actions 在 **60 天内没有任何活动**，相关的**定时 workflow** 会被自动禁用。

这意味着我每过两个月左右就需要手动去 GitHub 仓库页面重新启用这些 workflow。如果不这样做，channel 的自动更新就会中断。所以，目前最省心的方式还是在手机上设个提醒，每 60 天去重新开启一次。虽然不优雅，但胜在简单可靠。

## channel 的使用情况

起初我以为这个 channel 只有我自己在用，后来发现其实也有其他人通过它安装软件包。

其中被下载最多的包是 **opencode**。不过现在 opencode 已经在 GitHub Release 上提供了直接下载的渠道，以后通过 channel 安装的用户应该会越来越少。这也意味着我的 channel 维护压力会逐渐减轻——至少从 opencode 这个包的角度来看是这样。
