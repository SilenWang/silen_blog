---
title: 使用 DevPod 设置自部署的 Codespace
categories: Coding
date: 2025-10-28 05:49:54
tags: [DevPod, DevOps, 开发环境, Docker, Podman]
---

DevPod 是一个开源的开发环境管理工具，可以让你在任何 Kubernetes 集群或 Docker 主机上创建类似 GitHub Codespaces 的开发环境。本文将介绍如何使用 DevPod CLI 创建工作区，并详细解析 DevContainer 配置文件的编写。

<!-- more -->

## DevPod 介绍

DevPod 是由 Loft Labs 开发的开源工具，它允许开发者在任何基础设施上创建可重复、一次性的开发环境。与 GitHub Codespaces 类似，但 DevPod 是可以自托管的，可以在任何在本地机器、支持的云服务器商或 Kubernetes 集群上运行。

主要特性包括：
- 支持多种后端（Docker、Kubernetes、AWS EC2 等）
- 基于 DevContainer 标准
- 开发环境即代码
- 快速启动和销毁环境

## DevPod Provider 和工作区创建

### 安装 DevPod

Devpod 提供桌面版本，但是我现在是要在 FydeOS 上运行，所以安装 DevPod CLI。比较友好的是 Devpod CLI 几乎没有外部依赖，甚至SSH都有内建的，因此使用官方的命令就能安装到 FydeOS 的原生命令行。

```bash
curl -fsSL https://raw.githubusercontent.com/loft-sh/devpod/master/scripts/install.sh | sh
```

### 配置 Provider

Provider 是 DevPod 的后端驱动，定义了开发环境的运行位置。以 SSH 为例：

```bash
devpod provider add ssh --name amd -o HOST=AMD
```

上面的添加方式基于已有的 ssh 配置，对应 `~/.ssh/config` 中的内容是：

```
Host AMD
  HostName 192.168.0.2
  Port 22
  User user
  IdentityFile /home/user/.ssh/AMD_Key
```

另须注意，DevPod 是基于容器技术的工具，因此 SSH 的目标机器需要已经安装好 docker 或者 podman。

另外，使用 Ubuntu 的话，还要额外注意 docker 软件的来源，因为 Ubuntu 下可以通过 Snap 安装 docker，这个版本的 docker 因为 Snap 的限制，是不能访问用户 Home 目录下的隐藏文件夹的（`.`开头的那些），导致 Devpod 调用 docker 执行任何构建命令都会失败。

### 创建工作区

指定 provider 和项目路径即可启动工作区，Devpod会根据配置文件创建必要的工作环境，更重要的是，Devpod 会自动在工作区安装 web 版的 openvscode，启动该服务后，将服务的端口自动转发到本地，然后就可以快乐的写代码了，真正实现了项目介绍中说的，开发环境即配置的目标。同时，这些配置只需要服务器开 SSH 端口即可使用，可以说，真的相当方便。

```bash
devpod up --provider amd --source git https://github.com/your-username/your-repo
```

这将基于仓库中的 `.devcontainer` 文件夹下的配置文件创建开发环境。

## DevContainer 格式介绍

DevContainer 是开发环境配置的标准格式，本次示例中包含两个主要文件：`devcontainer.json` 和 `Dockerfile`。

### Dockerfile 配置

Dockerfile 其实就是 docker build 时使用的配置文件，DevPod会利用这个文件来创建工作区的容器。

需要注意的是，配置构建容器最好基于 mcr.microsoft.com/devcontainers 下的容器来构建，因为这些容器是根据一定的要求构建而成的，DevPod 的一些功能依赖容器内的文件和配置，如果使用其他容器，启动工作区会需要自行处理一些兼容性问题。

```dockerfile
# 使用 ms 官方基础镜像，含有vscode用户，避免devpod启动时有问题
FROM mcr.microsoft.com/devcontainers/base:bookworm

# 为vscode安装pixi
USER vscode

# 安装pixi，因为我的项目都用它管理依赖
RUN curl -fsSL https://pixi.sh/install.sh | sh

# 配置 git 用户信息，同时向环境注入我需要的alias内容
RUN git config --global user.name "Sylens" && \
  git config --global user.email "qiumin14@163.com" && \
  echo 'alias aider="aider --no-check-update --no-show-model-warnings --yes --no-auto-commits --model deepseek/deepseek-reasoner"' >> ~/.bashrc

ENV PATH="/home/vscode/.pixi/bin:${PATH}"
```

另外还需要注意的一点是，DevPod的运行逻辑是先构建容器，然后再根据 `devcontainer.json` 进行项目代码挂载和进一步配置，所以一切依赖项目代码的设置都是无法在容器构建阶段进行的。

### devcontainer.json 配置

`devcontainer.json` 定义了开发容器的元数据和 IDE 配置，我这里是个最基本的配置，只设置了一些环境变量以及最低限度的插件情况。

```json
{
  "name": "silen_blog",
  "build": { "dockerfile": "Dockerfile" },
  "postCreateCommand": "pixi run deploy",
  "containerEnv": {
    "AIDER_DARK_MODE": "${localEnv:AIDER_DARK_MODE}",
    "AIDER_CODE_THEME": "${localEnv:AIDER_CODE_THEME}",
    "DEEPSEEK_API_KEY": "${localEnv:DEEPSEEK_API_KEY}"
  },
  "customizations": {
    "vscode": {
      "settings": {
        "workbench.colorTheme": "Solarized Dark"
      },
      "extensions": [
        "naumovs.color-highlight",
        "ms-ceintl.vscode-language-pack-zh-hans",
        "tamasfe.even-better-toml"
      ]
    }
  }
}
```

须注意，`postCreateCommand` 只能写一项，如果需要执行的内容过多，建议写成shell脚本，放在 `.devcontainer` 目录然后用 `postCreateCommand` 调用。

## 后记

早在去年写 [ReviewGPT](https://github.com/SilenWang/ReviewGPT) 的时候（{% 在chatGPT的指导下写一个调用chatGPT的API进行文献解析的工具 ReviewGPT [在chatGPT的指导下写一个调用chatGPT的API进行文献解析的工具] %}），我就体会到，想法人人有，做出来最重要。

我本来是在改写ChromeOS下的终端App，希望为它加入端口转发，和在目标机器上自动安装 vscode web 并自动转发端口到本地的功能。然后做着做着，发现google 早就废弃 pNaCl 了（如此一来原生的 App 并不能访问本地端口了），于是想用 Golang 开发一个满足我要求的引用，继续做着做着，发现了 DevPod，除了它是基于容器的... 已经非常接近我的需求了... 害，那我还写啥呢... 先学 Devpod 怎么用吧...