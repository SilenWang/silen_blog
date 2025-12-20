---
title: git中的subtree与submodule
categories: Other
date: 2025-12-06 01:49:33
tags:
  - git
  - 版本控制
  - subtree
  - submodule
---

Emmmmm，今年官网供应商的合同到期了，于是... 多了一个项目要管... 这个项目依旧是前后端分离的，不同的是，官网有英文版本，且英文版本是前端项目的一个分支。我以前用过Submodule，以将同事写的独立模块整合到主项目中，但是这次我在跟AI请教后，选择了另外一种方式，Subtree

<!-- more -->

## 快速比较：Submodule vs Subtree

简单来说，Submodule 是在主仓库中创建一个指向子仓库特定提交的“链接”，而 Subtree 则是将子仓库的代码完全复制到主仓库的一个子目录中。两者各有优劣，下表可以帮助你快速了解它们的核心区别：

| 特性 | Submodule | Subtree |
|------|-----------|---------|
| 代码存放方式 | 主仓库仅保存子仓库的引用（commit hash） | 子仓库的代码完全复制到主仓库中 |
| 克隆后是否需要额外初始化 | 需要 `git submodule init && git submodule update` | 无需额外操作，所有文件已就绪 |
| 在主项目中直接修改子项目代码 | 必须进入子模块目录单独提交，步骤繁琐 | 可直接在主仓库目录中修改并提交 |
| 历史记录 | 子仓库历史与主仓库完全分离 | 子仓库的历史可选择性地合并（`--squash`）或完整保留 |
| 适合场景 | 子项目相对稳定，更新频率低；多个主项目共享同一个子项目 | 需要频繁在主项目中修改子项目；需引入同一仓库的多个不同分支 |

简而言之，`Subtree`实际上是把子项目代码完整克隆了一份，修改后提交的时候是将主项目里产生的修改记录同步一份到子项目，而子项目被修改时，也可以将这些修改同步回主项目，而`submodule`实际上是软链接，主项目中只是指定用哪个版本的子项目代码，提交代码回子项目会相对麻烦一些。

对我而言，一方面，我想要**将同一个仓库的两个不同分支（`master` 与 `en`）分别作为独立的子目录引入**；另一方面，我实际需要频繁修改子项目的代码，主项目实际是对这前端后端的整合，根据AI的建议，**Subtree**是操作更简单的组合。

## Subtree 的实际使用方式

下面是我在主项目中添加各个子项目的具体命令。注意，我使用了 `--squash` 参数，它可以将子项目的历史提交合并为一个，避免主项目历史被过多的子项目提交淹没。

### 1. 添加远程仓库
```bash
git remote add Backend https://github.com/user/Backend
git remote add Frontend https://github.com/user/Frontend
```

### 2. 将子项目以 subtree 形式引入
```bash
# 中文网站（master 分支）
git subtree add --prefix=app/Backend Backend master --squash
# 英文站独立分支（en 分支）
git subtree add --prefix=app/Backend_EN Backend en --squash

# 后端管理项目
git subtree add --prefix=app/Backend Backend_Admin aiming_med --squash
```

### 3. 日常同步操作
- **拉取子项目更新**：`git subtree pull --prefix=<目录> <远程名> <分支> --squash`
- **推送主项目内对子项目的修改**：`git subtree push --prefix=<目录> <远程名> <分支>`

当然为了简化操作，可以在Pixi中配置响应的命令

### 4. 注意事项
- 使用上述方式设置 subtree 后，子项目代码在本项目已有一份完整拷贝，因此，再次部署时，克隆主项目后**无需**再次运行 subtree add。
- 如果是为了开发而克隆，仍需添加远程仓库（如上第一步），以便后续推送更新信息到子项目。

## 后记

实际使用了Subtree两周，除了踩到个Devcontainer的git没有subtree模块的坑，其他方面倒是运行良好... 希望不要有什么新问题吧...