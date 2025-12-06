---
title: pixi也可以在github_action中使用了
categories: Coding
date: 2025-12-06 01:43:09
tags: ['github action', 'pixi']
---

在 CI/CD 流程中，依赖管理往往是决定构建效率与可靠性的关键因素。最近，我在一个静态网站部署流水线中尝试了 [setup-pixi](https://github.com/prefix-dev/setup-pixi) 这个 GitHub Action，

<!-- more -->


## 为什么选择 Pixi？

我最开始选择Pixi的时候，主要是出于因为它在使用 `conda` 的预编译资源基础上，具有比mamba更快的依赖计算速度，且能跟 `nodejs` 一样有依赖锁机制，适合进行快速的环境迁移部署，这些在软件经常变动的生物信息分析中都是非常需要的特性。

而在使用了一年多后，我不仅使用它管理生物信息分析的环境，也进一步用他来建立和管理开发环境，目前除了MySQL这种必须要系统常驻服务的没有办法处理外，Pixi完全可以安装所有我工作中需要用到的工具链。

使用单一的工具就可以进行各种工具链的管理，并且也内置简单的工作流设置，使用体验相当之好。

唯一的美中不足，就是之前根据领导要求使用Github Action做自动化的时候有那么一点点不便，因为 `pixi` 毕竟还是涉及虚拟环境，安装和运行的时候都需要进行一些额外的环境变量操作才能使用。

但是其实官方早就开发了`setup-pixi` [action](https://github.com/marketplace/actions/setup-pixi)，只是我一致不知道而已...

## 在 GitHub Actions 中使用 `setup-pixi`

`setup-pixi` 使用起来很简单，下面是一个例子：

```yaml
name: DEPLOY

on:
  push:
    branches:
      - main
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

      - name: Run Pixi Task
        run: |
          pixi run build
```

### 步骤解析

1. `Checkout repository`：使用 `actions/checkout` 获取最新代码，需要在 `Setup Pixi` 之前，毕竟 `pixi.toml` 和 `pixi.lock` 在项目里；
2. `Setup Pixi`：`prefix-dev/setup-pixi` Action 会下载并安装指定版本的 Pixi，并自动激活我们在`matrix` 中指定，并在 `pixi.toml` 中定义过的环境（本例中的 `website` 环境）；
3. `Run Pixi Task`：有了环境，就可以任意调用项目中已经定义过的task了。

整个过程简洁清晰，Pixi 环境保证了构建阶段依赖的一致性，一切交给它，定义一次，就不需要再做其他的设置了。


## 小结

通过在 GitHub Actions 中利用 `setup-pixi`，我们可以轻松地使用 `pixi` 的强大依赖管理能力，简化配置过程，减少人工设置的时间，并保证了从开发到生产环境的高度一致。