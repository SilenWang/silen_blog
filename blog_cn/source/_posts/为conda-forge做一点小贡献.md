---
title: 为conda_forge做一点小贡献
categories: Others
date: 2026-01-09 10:03:56
tags: ['conda-forge', 'tree_sitter_languages', 'PR', 'conda', 'feedstock']
---

Conda 上的终于有 `aider-chat` 的 conda包了，这样理论上可以通过 `pixi global` 来全局安装它了。但是实际安装会发现，其依赖项之一 `tree_sitter_languages` 却没有对应的 aarch64 版本，安装过程会因此失败。于是我就想到，我能不能靠AI来解决这个。

<!-- more -->

## 什么是 conda‑forge？

说来惭愧，在这之前，我一致不明白 conda 的 anaconda 和 conda‑forge 频道有什么区别... 这次翻文档才知道[Conda‑forge](https://conda-forge.org/) 是一个社区驱动的 Conda 软件包分发频道，而 Anaconda 是 conda 的官方资源（以及如果在商业公司随便使用Anaconda，可能会被警告和要求付费...）。
而 Conda‑forge 上的每一个软件包都对应一个 **feedstock** 仓库，其中包含了构建该包所需的全部文件，其中最关键的就是 **recipe**。

## 什么是 recipe？

Recipe 是指导 Conda 如何构建一个软件包的配置文件，很形象的被叫成了“菜谱”，具体可以见[官方文档](https://conda-forge.org/docs/maintainer/adding_pkgs)。

由于我本次并不是从头构建一个conda包，所以实际上这次并没有涉及修改 Recipe。

## 准备 feedstock 仓库

如前所述，每一个软件包都对应一个 **feedstock** 仓库，[`tree_sitter_languages‑feedstock的仓库`](https://github.com/conda-forge/tree_sitter_languages-feedstock) 就在conda-forge组织下。根据官方文档的建议，对仓库做贡献建议先Fork仓库到自己的账户，然后开分支修改，推送到原仓库。


## 修改文件

根据 Aider 对这个 feedstock 的解析，这个库并没有什么特别不支持的 aarch64 的地方，可能只是单纯建立的时候就没有想支持 aarch64，并且默认状况下，feedstock 也不会自动支持 aarch64 （大概是，用户太少没人权吧...）。因此，修改很简单，只需要改 `conda-forge.yml`，加上：

```yaml
provider:
  linux_aarch64: azure
```

之后，根据文档，我需要手动运行`conda-smithy render -c auto`，conda 会自动对构建脚本进行更新（真自动化，这就是专业的CI和CD么...）

## 向原项目提交 PR 以进行更新

最后就是将代码提交到自己仓库的fork，然后从github向原项目提交PR，PR会触发原项目库的CI，这个流程会自动完成构建和测试，测试通过，包的Maintainer接受PR后，包就能在conda上搜到了。

有了这个包，在Fydetab Duo上就能直接`pixi global install aider-chat`了，目前使用没有出现啥问题~