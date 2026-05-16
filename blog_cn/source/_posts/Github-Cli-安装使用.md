---
title: Ubuntu安装并使用Github CLI
categories: Tools
date: 2026-05-10 12:00:00
tags: ['GitHub CLI', 'gh', 'Ubuntu', '认证', '命令行', 'git', 'issue']
---

原来拿到新电脑总要配置下github的ssh key，来方便克隆和推送代码，但是又不经常去保存记录key，导致每次都重新创建，还是用multica才在issue区知道，github是除了github-cli的，可以快速进行认证...

<!-- more -->

## 安装方式

### 方式一：apt 安装（推荐）

Debian 13 和新版本的 Ubuntu的软件源已经收录了 gh，直接安装就行。

```bash
sudo apt update
sudo apt install gh
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