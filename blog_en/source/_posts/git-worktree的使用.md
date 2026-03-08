---
title: Using git worktree
categories: Coding
date: 2026-03-07 23:29:26
keywords: [git, git worktree, multiple working directories, development environment, isolation]
---

<!-- 摘要部分 -->
git worktree is a powerful feature provided by Git that allows you to check out multiple branches in the same repository simultaneously. Like cherry-pick, I hadn't heard of this feature before, but in today's era of extensive AI Agent applications, this feature will likely become as fundamental and essential as commit.

<!-- more -->

## What is git worktree

git worktree allows you to have multiple working directories in the same Git repository, with each working directory checked out to a different branch. This feature is extremely useful when you need to work on multiple branches simultaneously without frequent branch switching. However, for individual developers or small teams, one person typically doesn't have the mental capacity to work on multiple unrelated tasks at once, so this feature isn't used that often in practice...

But Vibe Coding is different - since you're not writing code yourself anyway, and to reduce potential merge conflicts, you would actually prefer to work on unrelated tasks in parallel where conflicts are less likely. This makes this feature a priority to apply.

## Basic Usage

### Adding a worktree

```bash
# Create a new worktree at the specified path and automatically switch to the target branch
# Example: Create a worktree for the feature branch in ../my-feature directory
git worktree add ../my-feature feature-branch
```

### Listing all worktrees

```bash
git worktree list
```

### Removing a worktree

```bash
# First remove the worktree
git worktree remove /path/to/worktree

# If there are uncommitted changes, you need to force removal
git worktree remove --force /path/to/worktree

# Prune invalid worktree records
git worktree prune
```

## Practical Use Cases

### Scenario 1: Parallel Development of Multiple Features

Suppose you're working on Feature A, but also need to fix a bug in Feature B:

```bash
# Main repository continues development on feature-A branch
git checkout feature-A

# Create a new worktree on bugfix branch to fix the bug
git worktree add ../bugfix-branch bugfix-branch

# Both directories can be open simultaneously without interfering with each other
```

### Scenario 2: Maintaining Release Versions

```bash
# Main repository develops new features
git checkout develop

# Create a separate worktree for the release version
git worktree add ../release-v1.0 v1.0-release
```

## Comparison with Similar Features

As you can see, the ultimate goal of worktree is to create multiple forks, then modify them separately to reduce interference while making subsequent merging easier.

To achieve these goals, there are other options available, each with its own advantages and disadvantages:

### 1. git worktree

**Advantages:**
- High disk space efficiency (shares .git directory)
- No recompilation needed when switching branches (important for large C++/Go projects)
- Fast creation/deletion
- Completely native, no additional tools required

**Disadvantages:**
- Operates within the same repository, need to pay attention to branch status
- Not suitable for complete isolation scenarios - only the git-managed code is isolated

### 2. Multiple Repository Clones

**Advantages:**
- Complete isolation, no interference
- Can work on different repository versions simultaneously
- Simple operation, low learning curve

**Disadvantages:**
- Wastes disk space (each clone has a complete .git)
- Needs to reinstall dependencies for each clone
- Cloning large projects takes time

**Use Case:** When you need completely isolated code versions, or to share the same version with team members.

### 3. devcontainer

**Advantages:**
- Completely isolated development environment (code, runtime, database)
- Can configure unified development toolchain
- Reproducible development environment

**Disadvantages:**
- Requires Docker runtime environment
- Slow first-time startup
- Configuration has a learning curve
- Overly heavyweight for simple scenarios

**Use Case:** When teams need unified development environments, or projects depend on complex toolchains.

## Comparison Summary

| Feature | git worktree | Multiple Clones | devcontainer |
|---------|-------------|-----------------|--------------|
| Disk Usage | Low (shared .git) | High | Medium |
| Isolation Level | Medium | High | High |
| Startup Speed | Fast | Fast | Slow |
| Configuration Complexity | None | None | High |
| Use Case | Multi-branch parallel development | Version isolation | Need to isolate more than just code |

## Closing Thoughts

During the writing of this blog post, cc has already embedded worktree as a built-in feature... However, as a understanding of the principles and as a reference for technical choices in other development scenarios, this blog should still have some reference value, right?....
