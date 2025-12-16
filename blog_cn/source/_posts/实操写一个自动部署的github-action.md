---
title: 实操写一个自动部署的github_action
categories: Coding
date: 2025-12-17 01:44:13
tags: ['github action', 'CD', '自动部署']
---

随着手上要维护的内容越来越多（3个官网，2个系统，1个小程序，全部是前后端分离且数据库独立的项目），出现了很多做1次不费事，但是穿插着做很多次非常乱的工作。我之前已经尝试过用Github进行 CI，这次学习和实践了一下 CD。

<!-- more -->

## 个人维护多项目的痛点

作为被迫单打独斗的打工人，我需要同时维护多个技术栈各异的项目。每个项目都有独立的代码库、构建流程和部署环境。以往的做法是：本地修改代码 → 手动运行构建脚本 → 通过 SCP 将产物上传到服务器 → 登录服务器执行替换和重启操作。这种流程重复性极高，且容易因操作遗漏或手误导致线上问题。

尤其是当修复一个紧急 bug 后，需要在多个项目间来回切换时，手动部署的耗时和心智负担会显著增加。因此，我希望将这套重复性的操作交给机器自动完成，从而让自己更聚焦于代码逻辑本身。

## 为什么需要 CD（持续部署）

CI（持续集成）负责代码的自动化构建与测试，确保每次提交都是可运行的；而 CD（持续部署）则在此基础上，将通过测试的代码自动部署到生产或预览环境。对于需要频繁更新、多项目并行的团队，CD 能极大减少重复的手动操作，降低人为错误，让开发者更专注于功能开发。

当然，我不是这样的团队，我只有一个人，但如前所述，一个可用的 CD 工作流依然能让我免于进行机械的重复工作，将精力放到代码问题本身的修复上。

## 编写一个 GitHub Actions CD 工作流

下面是我为维护的官网项目设计的部署工作流。它会在代码推送到 `main` 分支、且变更发生在 `app/Official_Site/**` 或 `app/Official_Site_EN/**` 目录时自动触发，同时也支持手动在 GitHub 界面点击运行（`workflow_dispatch`）。

```yaml
name: DEPLOY

on:
  push:
    branches:
      - main
    paths:
      # 只有此目录发生变化才触发
      - 'app/Official_Site/**' 
      - 'app/Official_Site_EN/**'
  workflow_dispatch:

jobs:
  build:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        environment: [website]

    steps:
      - name: Checkout repository
        uses: actions/checkout@v3
        with:
          ref: main

      - name: Setup Pixi
        uses: prefix-dev/setup-pixi@v0.9.3
        with:
          pixi-version: v0.59.0
          environments: ${{ matrix.environment }}

      - name: Build Website for ZH/EN
        run: |
          pixi run web_build
          mv app/Official_Site/dist front
          mv app/Official_Site_EN/dist front-en

      - name: Copy files ZH
        uses: appleboy/scp-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          source: "front"
          target: /opt/project/uploadPath

      - name: Copy files EN
        uses: appleboy/scp-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          source: "front-en"
          target: /opt/project/uploadPath

      - name: Deploy
        uses: appleboy/ssh-action@v1
        with:
          host: ${{ vars.PROD_IP }}
          username: root
          password: ${{ secrets.PROD_PASSWD }}
          port: 22
          script: |
            cd /opt/project
            rm -rf front.bak front-en.bak
            mv front front.bak
            mv front-en front-en.bak
            mv uploadPath/front front
            mv uploadPath/front-en front-en
            systemctl restart nginx
            
```

### 关键步骤说明

1. **触发条件（`on`）**
   - `push` 到 `main` 分支且仅当指定目录有变动，避免无关提交触发构建。
   - `workflow_dispatch` 提供手动触发入口，方便在需要时立即部署。

2. **代码检出（Checkout）**
   - 使用官方 `actions/checkout` 动作拉取代码，并指定 `ref: main` 确保构建基于最新提交。

3. **依赖安装（构建产物（Build Website））**
   - 采用 [Pixi](https://pixi.sh/) 作为跨语言包管理工具，这里用它内置的命令同时完成环境设置和静态网站构建

4. **文件传输（SCP）**
   - 使用 `appleboy/scp-action` 将构建产物分别上传到服务器的临时目录 `/opt/project/uploadPath`。
   - 服务器 IP 通过 `vars.PROD_IP`（项目变量）传递，密码则存放在 `secrets.PROD_PASSWD`（仓库机密）中，避免敏感信息暴露在代码里。

5. **部署与回滚（SSH）**
   - 通过 `appleboy/ssh-action` 登录服务器执行替换操作。
   - 关键逻辑：先将当前运行的 `front` 和 `front-en` 目录备份为 `.bak` 后缀，再从上传目录移入新版本，最后重启 Nginx。
   - 备份操作提供了快速回滚的能力——只需将备份目录移回即可恢复上一版本。

## 总结

通过这个简单的 GitHub Actions 工作流，我成功将静态官网的部署时间从每次约 10+ 分钟的手动操作压缩到 2‑3 分钟的自动执行，并且完全杜绝了因操作遗漏导致的问题。

CD 并不是大团队的专利，个人项目同样可以通过自动化获得显著的效率提升。如果你也在维护多个需要频繁更新的项目，不妨花一小时配置一套属于自己的持续部署流水线，相信它会给你的开发体验带来质的飞跃。
