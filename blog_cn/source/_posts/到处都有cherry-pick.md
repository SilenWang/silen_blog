---
title: 到处都有 cherry-pick
categories: Coding
date: 2026-01-20 00:05:43
tags:
  - Git
  - 版本控制
  - "cherry-pick"
---

之前想给 Fydetab Duo Wiki 做点贡献，但是为了本地编译博客预览，我需要向项目中加入 pixi 或者其他配置，这些内容是不适合提交到原项目的。所以学习了一下如何选择性提交，也就是 `git cherry-pick`。

<!-- more -->

## Git 的 cherry-pick

在 Git 版本控制系统中，`cherry-pick` 命令用于将某个特定的提交（commit）应用到当前分支，而不需要合并整个分支。这很适合我现在的需求：我 Fork 了源代码，创建了一个用于修改的分支（dev），分支的前两个 commit 是我进行本地配置用的文件，后面的才是对文档内容修改的部分。当我完成修改准备提交时，我可以基于未修改的原分支，再创建一个提交专用分支，然后利用 `cherry-pick` 来选择要纳入的后三个提交：

```bash
# 将指定提交应用到当前分支
git cherry-pick COMMIT3 COMMIT4 COMMIT5

# 如果遇到冲突，解决后继续
git cherry-pick --continue

# 或者放弃 cherry-pick
git cherry-pick --abort
```

然后再使用这个分支提交 PR，就不会把不必要的内容提交到原项目了。

同时，上游合并了我的 PR 后，我也可以再从上游拉取更改到 dev 分支，这样就能继续更新，然后再修改和使用合并专用分支。

## 在数据分析中的 cherry-picking

当时看到 `cherry-pick` 的时候觉得还挺有趣的，因为这个词之前在数据分析相关的讨论中也见过，其意思比较负面，即"选择性展示"。说实在的，这种事情工作后干的还真不少... 也是一种统计学的"魅力"时刻了... 

另外，问了一下 AI，cherry pick 的说法似乎由来已久... Emmmmmm... 果然是太阳底下无新鲜事...