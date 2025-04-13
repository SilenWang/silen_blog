---
title: Github Actions的基本使用
categories: Script
date: 2025-04-13 16:40:47
tags: ['github', 'github actions']
---

在现代软件开发中，持续集成（CI）和持续部署（CD）是提高开发效率和质量的关键实践。[GitHub Actions](https://github.com/features/actions) 作为 GitHub 提供的强大自动化工具，能够帮助开发者轻松实现 CI/CD 以及更广阔的自动化工作流程。

<!-- 摘要部分 -->
<!-- more -->

## GitHub Actions 的作用
GitHub Actions 是一个强大的自动化平台，允许开发者在 GitHub 仓库中定义和运行工作流。它可以根据特定的事件（如代码推送、拉取请求、标签发布等）自动执行一系列任务，例如代码构建、测试、部署等。通过 GitHub Actions，开发者可以实现从代码提交到部署的全流程自动化，从而节省时间和精力，减少人为错误。


## GitHub Actions 的工作流结构
GitHub Actions 的工作流由多个部分组成，其中最重要的部分是 jobs 和 steps。

### Jobs
一个工作流可以包含多个 jobs，每个 jobs 是一个独立的执行单元，可以在不同的环境中运行。

```yaml
jobs:
  build: # 工作名称
    runs-on: [ubuntu-latest] # 运行环境，这里指定的是github hosted runners
```

runs-on 指定了运行该 jobs 的环境，这里我们使用github准备的最新ubuntu环境。

### Steps
每个 jobs 可以包含多个 steps，每个 steps 是一个独立的任务，可以执行一系列命令或调用预定义的动作（actions）。预定义的动作可以在[github marketplace](https://github.com/marketplace)搜索。

```yaml
- name: Get current tag # 步骤名称
  id: get_tag # 步骤编号
  uses: zingimmick/github-action-get-current-tag@v1 # 使用特定的已编写好的步骤
```

## 变量传递

就如同我过去写的生物信息分析流程，不同的步骤之间，常常需要进行信息传递才能完成一项工作。yaml作为简单的配置文件，自然无法实现变量传递，因此我们需要以来github action框架中设置的一些环境变量来完成信息传递。


### 跨Steps信息传递

在steps中执行代码时，可以通过向`$GITHUB_OUTPUT`这个特殊环境变量追加内容的方式来实现steps间信息传递的作用：
```yaml
- name: Save info
    id: save_info
    run: |
        echo "tag=3" >> $GITHUB_OUTPUT
- name: Use info
    id: use_info
    run: |
        echo ${{ steps.save_info.outputs.tag }}
```

### 跨Jobs信息传递

在jobs中，可以规定将特定步骤的输出传递到jobs输入，方法如下：

```yaml
jobs:
    get_info:
        steps:
            - name: Save info
                id: save_info
                run: |
                    echo "tag=3" >> $GITHUB_OUTPUT
    outputs:
      tag: ${{ steps.save_info.outputs.tag }}

    use_info:
        needs: [get_info] # 设定依赖，以获得其他jobs的信息
        steps:
            - name: Save info
                id: save_info
                run: |
                    echo ${{ needs.get_info.outputs.tag }}
```


## 编译项目并进行发布的例子

```yaml
name: BUILD

on:
  push:
    tags:
      - "v*" # 收到版本tag触发

jobs:
  build:
    runs-on: [ubuntu-latest]
    permissions: # 需要写release的权限
      contents: write

    steps: 

      - name: Get current tag
        id: get_tag
        uses: zingimmick/github-action-get-current-tag@v1
  
      - name: Checkout code
        uses: actions/checkout@v4
        with:
          fetch-depth: 0
          ref: ${{ steps.get_tag.outputs.tag }}

      - name: Build frontend
        run: pixi run front_prod
        
      - name: Build backend
        run: pixi run backend

      - name: Pack to zip
        run: pixi run release

      - name: Release
        uses: ncipollo/release-action@v1.15.0
        with:
          draft: false
          generateReleaseNotes: true  #自动生成发行说明。
          artifacts: '${{ github.workspace }}/release/*.zip'
          tag: ${{ steps.get_tag.outputs.tag }}
          name: Release ${{ steps.get_tag.outputs.tag }}
```