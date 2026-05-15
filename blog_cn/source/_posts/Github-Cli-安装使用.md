---
title: Ubuntu安装并使用Github CLI
categories: Tools
date: 2026-05-10 12:00:00
tags: ['GitHub CLI', 'gh', 'Ubuntu', '认证', '命令行', 'git', 'issue']
---

[GitHub CLI](https://cli.github.com/)（简称 `gh`）是 GitHub 官方推出的命令行工具，让你无需离开终端就能完成 GitHub 上的大部分操作。更棒的是，一旦完成认证，`git` 命令也会自动复用这份凭证，无需额外配置 SSH。

<!-- 摘要部分 -->
<!-- more -->

## 安装方式

### 方式一：apt 安装（推荐）

GitHub CLI 官方提供了 apt 仓库，可以直接用 `apt` 安装：

```bash
# 安装必要的依赖
sudo apt install curl

# 添加 GitHub CLI 的 GPG 密钥和仓库
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null

# 更新并安装
sudo apt update
sudo apt install gh
```

验证安装：
```bash
gh --version
```

### 方式二：snap 安装

如果你偏好 snap 包，也可以一行命令完成安装：

```bash
sudo snap install gh
```

## 登录认证

安装完成后，第一步就是登录 GitHub 账号：

```bash
gh auth login
```

执行后会进入交互式流程：

1. 选择登录方式：选择 **GitHub.com**
2. 选择协议：选择 **HTTPS**
3. 选择认证方式：选择 **Login with a web browser**，然后复制出现的一次性验证码
4. 在浏览器中粘贴验证码并授权
5. 终端确认后会显示认证成功

也可以使用 Token 进行非交互式登录：

```bash
gh auth login --with-token < my_token.txt
```

## 认证后 git 的联动效果

`gh auth login` 完成认证后，`git` 会自动利用 `gh` 的凭证管理器，这意味着：

- 对于基于 HTTPS 的远程仓库，`git push`、`git pull` 等操作不再需要每次都输入用户名和密码
- 无需手动配置 SSH key 就能进行代码推送
- 所有的认证信息通过 `gh` 统一管理

```bash
# 克隆私有仓库（无需额外配置）
git clone https://github.com/yourname/private-repo.git

# 推送代码（自动复用 gh 的凭证）
git push origin main
```

## GitHub CLI 的常用功能

### 操作 Issue

```bash
# 列出当前仓库的 Issue
gh issue list

# 查看 Issue 详情
gh issue view 123

# 创建新的 Issue
gh issue create --title "bug report" --body "describe the bug"

# 关闭 Issue
gh issue close 123
```

### 操作 Pull Request

```bash
# 列出 PR
gh pr list

# 创建 PR
gh pr create --title "fix bug" --body "fix description"

# 查看 PR 详情
gh pr view 123

# 合并 PR
gh pr merge 123
```

### 仓库操作

```bash
# 查看仓库信息
gh repo view

# 在浏览器中打开仓库
gh repo view --web

# Fork 一个仓库
gh repo fork
```

## 总结

GitHub CLI 让终端操作 GitHub 变得异常流畅。一次 `gh auth login` 搞定认证后，不仅 `gh` 本身的命令可以用了，`git` 命令也免去了反复输入凭证的麻烦。对于经常与 GitHub 打交道的开发者来说，这是非常值得安装的工具。
