---
title: 为conda_forge做一点小贡献
categories: Others
date: 2026-01-09 10:03:56
tags: ['conda-forge', 'tree_sitter_languages', 'PR', 'conda', 'feedstock']
---

在 conda‑forge 生态中，许多软件包默认只提供 x86_64 架构的预编译二进制文件。然而随着 ARM64（aarch64）设备的普及，为更多软件包添加对 ARM 平台的支持也成了社区的重要需求。本文将以 [`tree_sitter_languages‑feedstock`](https://github.com/SilenWang/tree_sitter_languages-feedstock) 为例，介绍如何为这类没有明确限制平台的包添加 aarch64 支持，并最终向原项目提交 Pull Request。

<!-- more -->

## 什么是 conda‑forge？

[Conda‑forge](https://conda-forge.org/) 是一个社区驱动的 Conda 软件包分发平台。它通过一套自动化的构建系统，将上游的源代码打包成能在 Conda 环境中直接安装的二进制文件。Conda‑forge 上的每一个软件包都对应一个 **feedstock** 仓库，其中包含了构建该包所需的全部文件，其中最关键的就是 **recipe**。

## 什么是 recipe？

Recipe 是指导 Conda 如何构建一个软件包的“菜谱”，它通常包含以下几个文件：

- `meta.yaml`：定义包的名称、版本、依赖关系、构建脚本等核心信息。
- `build.sh`（Linux/macOS）或 `bld.bat`（Windows）：执行实际的编译、安装命令。
- 可选的其他文件，如补丁、资源文件等。

Recipe 中可以通过 **跳过条件（skip）** 来控制包在哪些平台上构建。例如，如果某个软件目前只能在 64 位 Linux 上运行，recipe 中可能会包含：

```yaml
build:
  skip: True  # [not (linux and x86_64)]
```

这行代码表示“除非是 Linux x86_64 平台，否则跳过构建”。这种限制正是我们添加新架构时需要修改的地方。

## 检查现有 feedstock 的构建限制

以 [`tree_sitter_languages‑feedstock`](https://github.com/SilenWang/tree_sitter_languages-feedstock) 为例，假设它目前没有为 aarch64 提供预编译包。我们首先需要查看它的 recipe 中是否含有上述平台限制。

1. **Fork 并克隆 feedstock 仓库**

   ```bash
   git clone https://github.com/<你的用户名>/tree_sitter_languages-feedstock.git
   cd tree_sitter_languages-feedstock
   ```

2. **查看 `recipe/meta.yaml` 中的 `skip` 字段**

   打开 `recipe/meta.yaml`，寻找 `build:` 段落下的 `skip:` 行。如果该字段不存在，说明 recipe 本身没有显式限制平台，此时 aarch64 支持可能仅需在 CI 配置中开启。如果存在类似

   ```yaml
   skip: True  # [not (linux and x86_64)]
   ```

   的代码，则需要将其修改为同时允许 aarch64 架构：

   ```yaml
   skip: True  # [not (linux and (x86_64 or aarch64))]
   ```

   这里的 `# [ ]` 是 [Conda 构建变体（variant）](https://docs.conda.io/projects/conda-build/en/latest/resources/variants.html)的注释语法，它会在构建时根据目标平台自动求值。

## 启用 CI 对 aarch64 的构建

仅仅修改 `skip` 条件还不够，我们还需要确保 conda‑forge 的 CI 系统会为 aarch64 平台触发构建。这通常需要在 `.ci_support/` 目录下存在对应的配置文件，或者更简单的方式是让 conda‑forge 的自动迁移机制来生成它们。

对于社区贡献者，最直接的做法是 **向 feedstock 仓库提交一个包含以下更改的 PR**：

1. **修改 `conda-forge.yml`（如果存在）**  
   在 feedstock 根目录下的 `conda-forge.yml` 中添加（或确保已存在）以下内容：

   ```yaml
   provider:
     linux_aarch64: azure
   ```

   这表示使用 Azure Pipelines 来执行 Linux aarch64 的构建任务。

2. **创建或更新 CI 配置文件**  
   实际上，conda‑forge 的 CI 系统会在 PR 中自动生成 `.ci_support/` 下的配置文件。我们只需确保 PR 中包含了对 `recipe/meta.yaml` 的修改即可。一旦 PR 被打开，conda‑forge 的机器人（如 `@conda-forge-admin`）会自动重新渲染 feedstock，生成所有必要的 CI 文件。

## 实战步骤：以 tree_sitter_languages‑feedstock 为例

假设经过检查，该 feedstock 的 `recipe/meta.yaml` 中并没有 `skip` 字段，说明平台限制可能来自于其它地方（例如构建脚本中隐含的架构判断）。在这种情况下，我们依然可以尝试直接开启 aarch64 构建。

- **步骤 1：在本地创建新分支**

  ```bash
  git checkout -b add-aarch64-support
  ```

- **步骤 2：修改 `recipe/meta.yaml`（可选）**  
  如果不存在 `skip` 字段，我们也可以显式地加入一个允许 aarch64 的 skip 语句，这样可以避免未来被误限制。在 `build:` 段落下添加：

  ```yaml
  build:
    skip: True  # [not (linux and (x86_64 or aarch64))]
  ```

  如果 `build:` 段落已存在其他内容，请将 `skip:` 行放在其中。

- **步骤 3：提交并推送**

  ```bash
  git add recipe/meta.yaml
  git commit -m "Enable aarch64 builds"
  git push origin add-aarch64-support
  ```

- **步骤 4：在 GitHub 上创建 Pull Request**  
  访问你 fork 后的仓库页面，GitHub 通常会提示你刚刚推送的分支，点击 “Compare & pull request” 即可。在 PR 描述中简要说明此次更改的目的，例如：

  > This PR adds aarch64 (ARM64) support to the `tree_sitter_languages` package by adjusting the skip condition in `meta.yaml`. The change will allow the conda‑forge CI to build and distribute binaries for Linux aarch64 platforms.

## 等待 CI 与合并

提交 PR 后，conda‑forge 的 CI 流水线会自动开始构建，你可以在 PR 下方的检查列表中看到 `linux_aarch64` 的构建任务。如果构建成功，通常会有维护者来 review 并合并你的 PR。如果构建失败，则需要根据错误日志调整 recipe（例如可能需要在 `build.sh` 中添加针对 ARM 架构的编译选项）。

## 向原项目提交 PR 以进行更新

有时，你希望上游的源码仓库也注意到 ARM 平台的支持。例如，`tree_sitter_languages` 本身可能没有在 CI 中测试 ARM。此时你可以向上游仓库提交一个类似的 PR，帮助他们在源码层面更好地适配 ARM。

1. **找到上游仓库**  
   通常在 feedstock 的 `recipe/meta.yaml` 中会有 `source: url` 字段，指向源码地址。

2. **按照上游的贡献指南**，在源码中添加 ARM 相关的 CI 配置（如 GitHub Actions 中增加 `arm64` runner）或者在文档中注明对 ARM 的支持。

## 小结

为 conda‑forge 软件包添加 aarch64 支持并不复杂，核心步骤是：

1. **检查 recipe 中的平台限制**，必要时修改 `skip` 条件。
2. **确保 CI 配置包含 aarch64 提供器**（如 `conda-forge.yml`）。
3. **提交 PR 并等待 CI 构建成功**。

通过这样的贡献，你可以帮助更多用户在 ARM 设备上方便地使用 conda 生态中的工具，也让 conda‑forge 社区变得更包容、更强大。

---

撰写本文时参考了 [conda‑forge 官方文档](https://conda-forge.org/docs/maintainer/adding_pkgs.html) 以及多个已有的 ARM64 迁移 PR。如果你在操作中遇到问题，欢迎在 conda‑forge 的 [GitHub Discussions](https://github.com/conda-forge/conda-forge.github.io/discussions) 或 [Gitter 聊天室](https://gitter.im/conda-forge/conda-forge.github.io) 中寻求帮助。
