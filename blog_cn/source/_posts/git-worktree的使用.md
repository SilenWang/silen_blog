---
title: git worktree 的使用
categories: Coding
date: 2026-03-07 23:29:26
keywords: [git, git worktree, 多工作目录, 开发环境, 隔离]
---

<!-- 摘要部分 -->
git worktree 是 Git 提供的强大功能，允许在同一仓库中同时检出多个分支。这个功能和cherry-pick一样，我之前完全没听过，但是在变成Agent高度应用的今天，这个功能以后大概会变成跟commit一样的基础必学的功能。

<!-- more -->

## 什么是 git worktree

git worktree 允许你在同一个 Git 仓库中同时存在多个工作目录，每个工作目录可以 checkout 不同的分支。这个功能对于需要同时在多个分支上同时工作，但不想频繁切换的场景非常有用。但是，对于个人开发或者小团队来说，一个人其实不会脑容量高到可以同时开发几个不相干的内容，所以实际用到的时候不会那么多...

但是Vibe Coding就不一样了，反正都不是自己写，且为了降低合并代码时可能的冲突，反而会优先选择冲突可能小的不相干内容并行做，因此这个功能就变成要优先应用的东西了。

## 基本用法

### 添加工作目录

```bash
# 在指定路径创建新的 worktree，并自动切换到目标分支
# 示例：在 ../my-feature 目录创建 feature 分支的工作目录
git worktree add ../my-feature feature-branch
```

### 列出所有 worktree

```bash
git worktree list
```

### 删除工作目录

```bash
# 先移除 worktree
git worktree remove /path/to/worktree

# 如果有未提交的更改，需要强制移除
git worktree remove --force /path/to/worktree

# 清理无效的 worktree 记录
git worktree prune
```

## 实际使用场景

### 场景一：并行开发多个功能

假设你正在开发功能 A，同时需要修复功能 B 的 bug：

```bash
# 主仓库继续在 feature-A 分支开发
git checkout feature-A

# 新建 worktree 在 bugfix 分支修复 bug
git worktree add ../bugfix-branch bugfix-branch

# 两个目录可以同时打开，互不影响
```

### 场景二：保持发布版本

```bash
# 主仓库开发新功能
git checkout develop

# 为发布版本创建独立的 worktree
git worktree add ../release-v1.0 v1.0-release
```

## 其他类似功能对比

可以看出，worktree的最终目的是创建多个fork，然后在这些fork上分别修改，减少互相干扰，同时要方便后续的合并。

如果要实现这几点，其实还有其他的选择，不同的方案各有优劣：

### 1. git worktree

**优点：**
- 磁盘空间利用率高（共享 .git 目录）
- 切换分支无需重新编译（如大型 C++/Go 项目）
- 创建/删除速度快
- 完全原生，无需额外工具

**缺点：**
- 同一仓库内操作，需要注意分支状态
- 不适合完全隔离的场景，被隔离的只是被git管理的代码

### 2. 多次克隆仓库

**优点：**
- 完全隔离，互不影响
- 可以同时在不同仓库版本工作
- 操作简单，理解成本低

**缺点：**
- 磁盘空间浪费（每个克隆都有完整的 .git）
- 每次克隆都需要重新安装依赖
- 大型项目克隆耗时

**适用场景：** 需要完全隔离的不同版本代码，或与团队成员共享同一版本时。

### 3. devcontainer

**优点：**
- 完全隔离的开发环境（代码、运行环境、数据库）
- 可配置统一的开发工具链
- 可复现的开发环境

**缺点：**
- 需要 Docker 运行环境
- 首次启动较慢
- 配置有一定学习成本
- 对于简单场景过于重量级

**适用场景：** 团队协作需要统一开发环境，或项目依赖复杂的工具链。

## 对比总结

| 特性 | git worktree | 多次克隆 | devcontainer |
|------|-------------|----------|--------------|
| 磁盘占用 | 低（共享.git） | 高 | 中 |
| 隔离程度 | 中 | 高 | 高 |
| 启动速度 | 快 | 快 | 慢 |
| 配置复杂度 | 无 | 无 | 高 |
| 适用场景 | 多分支并行开发 | 版本隔离 | 代码以外的内容也需要隔离 |

## 后话

在写这个博客的过程中，cc已经将worktree作为内置功能嵌入了，不过... 作为原理了解和其他开发场景的技术选择参考，这个博客应该还是有那么点参考价值的，吧....
