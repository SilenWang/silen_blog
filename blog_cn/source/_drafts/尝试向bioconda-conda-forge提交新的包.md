---
title: 尝试向bioconda_conda_forge提交新的包
categories: Bioinformatics
date: 2026-01-11 17:15:43
tags: []
---

再上次尝试了对conda-forge的包做小贡献之后，我想继续来点更进阶的：尝试将 [singler‑py](https://github.com/SingleR-inc/singler-py) 这个 Python 包发布到 conda 生态中，结果不做不知道，一做... 还有点麻烦...

<!-- more -->

## 罗马并不能一天建成

我原本以为，设置好recipe，一次提交，通过审核就行了，但是发现了如下问题：

- 原本打算直接提交到 **conda‑forge**，却在阅读贡献指南时发现：专门面向生物信息学的软件应当优先提交到 **Bioconda** 频道。
- 之后阅读**Bioconda**的文档并尝试逐步尝试，发现 singler‑py 依赖的不少包都未进入任何 conda 频道，AI告诉我官方并不建议在同一个 PR 中提交多个新包，因此我只能从依赖树的最底层开始一个个来。
- 逐级往下解析，发现BiocPy下的一系列包都是没有进入Conda生态的... 梳理起来，从最底层的 biocutils开始，然后基础的 `biocframe`，再到数据结构的`summarizedexperiment`，最后才能到实际上应用层的`singler‑py`，如果每个提交都要2~3日审核... 实际上弄完至少是2个工作周吧...

嘛，反正都是学习，就一步步来好了...

## 提交到 singler‑py 位置需要的各步骤依赖

## conda-forge和Bioconda不同的recipe工作流

### conda‑forge 与 Bioconda 简介

### conda‑forge 提交工作流
[conda‑forge](https://conda-forge.org/) 是一个由社区主导的 conda 软件包仓库，它覆盖了绝大多数通用领域的开源软件。任何人都可以通过 GitHub 向它的 [staged‑recipes](https://github.com/conda-forge/staged-recipes) 仓库提交新的配方（recipe），经过自动化检查和维护者审核后，新包就会出现在 `conda‑forge` 频道中，供全球用户通过 `conda install -c conda-forge <package>` 安装。

### Bioconda 提交工作流
[Bioconda](https://bioconda.github.io/) 是一个专注于生物信息学软件的 conda 频道。它建立在 conda‑forge 的基础设施之上，所有从 Bioconda 安装的包都会自动依赖 conda‑forge 中的通用库（如 numpy、pandas 等）。Bioconda 的贡献流程与 conda‑forge 类似，但拥有自己独立的配方仓库 [bioconda‑recipes](https://github.com/bioconda/bioconda-recipes) 和一套专门为生物信息学软件设计的构建工具（`bioconda‑utils`）。

**重要原则**：如果一个软件的主要用途属于生物信息学范畴（例如单细胞分析、基因组组装、序列比对等），应当优先提交到 Bioconda；如果是通用工具（如文本处理、网络库、数学库等），则更适合提交到 conda‑forge。

### 2. 向 conda‑forge 提交新包的简要步骤

1. **Fork 仓库**  
   访问 [conda‑forge/staged‑recipes](https://github.com/conda-forge/staged-recipes) 并点击 “Fork” 按钮，将仓库复制到自己的 GitHub 账号下。

2. **准备本地环境**  
   克隆你 fork 后的仓库，并创建一个新分支：
   ```bash
   git clone https://github.com/<你的用户名>/staged-recipes.git
   cd staged-recipes
   git checkout -b add-<包名>
   ```

3. **编写配方文件**  
   在 `recipes/<包名>/` 目录下创建 `meta.yaml`（可参考其他现有配方）。该文件需要包含包的元信息（名称、版本、描述）、依赖列表、构建指令等。确保遵循 [conda‑forge 的文档规范](https://conda-forge.org/docs/maintainer/adding_pkgs.html)。

4. **本地验证**  
   使用 `conda smithy` 工具在本地运行 lint 和构建测试（需先安装 `conda‑smithy`）：
   ```bash
   conda install conda-smithy
   conda smithy recipe-lint recipes/<包名>/
   conda build recipes/<包名>/
   ```

5. **提交 PR**  
   将修改推送到你的 fork，然后在 GitHub 界面向 `conda‑forge/staged‑recipes` 发起 Pull Request。CI 会自动运行多平台构建测试，维护者会在评论中提出修改意见，直至所有检查通过后合并。

### 3. 向 Bioconda 提交新包的简要步骤

1. **Fork 仓库**  
   访问 [bioconda/bioconda‑recipes](https://github.com/bioconda/bioconda-recipes) 并 Fork 到自己的账号。

2. **准备本地环境**  
   克隆仓库并创建分支（同样建议分支名称包含包名）：
   ```bash
   git clone https://github.com/<你的用户名>/bioconda-recipes.git
   cd bioconda-recipes
   git checkout -b add-<包名>
   ```

3. **编写配方文件**  
   在 `recipes/<包名>/` 目录下创建 `meta.yaml`。Bioconda 的配方格式与 conda‑forge 基本相同，但有一些针对生物信息学软件的额外字段（例如 `extra: bioc‑install`）。可以参考仓库中已有的生物信息包配方。

4. **使用 bioconda‑utils 验证**  
   Bioconda 提供了专门的工具链来测试配方：
   ```bash
   # 安装 bioconda-utils（建议在独立 conda 环境中）
   conda create -n bioconda-utils bioconda-utils
   conda activate bioconda-utils

   # 运行 lint 检查
   bioconda-utils lint recipes/<包名>/

   # 在本地 Docker 容器中执行构建测试（可选，但推荐）
   bioconda-utils build --docker --packages <包名> .
   ```

5. **提交 PR**  
   推送到你的 fork 并向 `bioconda/bioconda‑recipes` 发起 PR。Bioconda 的 CI 会执行更严格的生物信息学软件兼容性测试，同样需要等待维护者审核并批准。

## 小发现

检查包构建文件的语法时，发现 conda现在也引入并行了，原来下载各个平台的包配置文件是串行的，大陆这即使用了加速能等个好几分钟，更不说后面的依赖计算了，现在引入了并行机制，还是比原来强了点（虽然我还是会继续用Pixi）

另外，印象中Pixi使用了更新的包构建工具，也许速度会比现在的conda生态更快，之后如果再提交别的，我一定试试。

## 后记

还在天津工作的时候（7年前了都...），就想过自己编译conda包，这样就能解决没有Root权限安装生信软件的问题了。只不过当时的文档是一个字也看不懂... 更别提写了。现在有了高度自动化的工具，也有了 AI 答疑... 确实可以来做点贡献了。