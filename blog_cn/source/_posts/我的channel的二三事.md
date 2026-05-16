---
title: 我的 channel 的二三事
date: 2026-05-08 07:00:00
tags:
  - channel
  - conda
  - GitHub Actions
categories:
  - 技术杂谈
---

## 前言

我在 [prefix.dev](https://prefix.dev) 上维护了一个 conda channel，用来分发一些自己打包的软件。随着时间推移，这个 channel 的维护工作量逐渐超出了我的预期，在这里记录一下遇到的一些问题。

<!-- more -->

## 配额限制

prefix.dev 上的 channel 并非无限的资源。每个 channel 的存储空间和包数量都有配额限制。当我往 channel 里上传的包越来越多时，很快就触及了配额上限。

目前我的做法是定期检查 channel 里的旧包，手动删除那些不再需要或者已经被新版本替代的包。这个过程完全是手工操作，比较繁琐，但暂时没有找到更好的自动化方案。

## GitHub Actions 的定时任务困境

我依赖 GitHub Actions 来自动构建和更新 channel 中的包。但 GitHub 有一个策略：如果一个仓库的 Actions 在 **60 天内没有任何活动**，相关的定时 workflow 会被自动禁用。

这意味着我每过两个月左右就需要手动去 GitHub 仓库页面重新启用这些 workflow。如果不这样做，channel 的自动更新就会中断。

我尝试过一些方案来规避这个问题：
- 添加一个 dummy commit 来触发 workflow？—— 这本质上还是手动操作。
- 用外部服务定期 ping 仓库？—— 增加了额外的依赖和复杂度。

目前最省心的方式还是在手机上设个提醒，每 60 天去重新开启一次。虽然不优雅，但胜在简单可靠。

## channel 的使用情况

起初我以为这个 channel 只有我自己在用，后来发现其实也有其他人通过它安装软件包。

其中被下载最多的包是 **opencode**。不过现在 opencode 已经在 GitHub Release 上提供了直接下载的渠道，以后通过 channel 安装的用户应该会越来越少。这也意味着我的 channel 维护压力会逐渐减轻——至少从 opencode 这个包的角度来看是这样。

## 后记

维护一个 conda channel 虽然有不少琐碎的麻烦事，但整体来说还是一个很有价值的经历。如果你也在考虑搭建自己的 channel，希望这篇小文能帮你提前了解可能会遇到的坑。

如果你也在维护 conda channel，欢迎交流经验！
