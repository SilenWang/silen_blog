---
title: 终于用 pixi 搞定了用于自动构建的 recipe 和 action：踩坑与解决方案
categories: Others
date: 2026-02-02 20:51:28
tags: ['pixi', 'recipe', 'github action', 'rattler', 'conda']
---

新工具总是这样：好用的地方让人惊艳，但尚未完善的部分却让人头疼。最近我尝试使用 pixi 和 rattler‑build 搭建一个自动化构建系统，用于定期打包并上传 opencode 到 prefix.dev。整个过程耗时约 6 小时，期间遇到了不少预料之外的问题。

<!-- more -->

## 背景

之前我已经成功为 opencode 创建了 recipe 并手动上传到了 prefix.dev。但 opencode 的更新频率相当高，我不可能每次都手动跟进。于是我想构建一个类似 bioconda 的自动化系统：存放 recipe，并通过 GitHub Action 自动检测新版本、构建并上传 conda 包。

感觉这工作应该不难？毕竟 recipe 已经有了，写个 Action 应该易如反掌，更何况还有 AI 辅助。然而这其中碰到的坑比我想的多的多的多...

## 主要挑战与解决方案

### 1. recipe 的 context 中不能使用 if‑then 判断

在 rattler‑build 的 recipe 中，`context` 部分用于定义变量，但这里**不能**像build里那样使用条件判断。AI 也是说 context 的部分只能静态的定义变量。然而这个说法... 对也不对。对的是，在 yaml 的语法内，确实是只能定义变量的，但是... rattler-build 其实是会先用 jinjia（是的又是jijia）来渲染出一个配置文件的，因此支持的 minijinjia 函数是可以实现判断的：

```yaml
context:
  architect: ${{ "linux‑x64" if platform == "linux‑64" else "linux‑arm64" }}
```

### 2. bump‑recipe 与 build 部分的兼容性问题

rattler‑build 的 `bump‑recipe` 功能可以自动更新 recipe 中的版本和 sha256，这非常方便，我可以不用在后续的 workflow 中用额外的 step 去进行一些操作了。但是，这个功能要求 recipe 的 `build` 部分不能包含条件判断。我之前的手动 recipe 在 build 部分使用了 if‑then，导致无法使用 bump‑recipe 功能。

于是我在更新后的 recipe 中全部使用 minijinjia 的函数， 不再使用 yaml 支持的特殊语法了。

### 3. AI 关于环境变量引用的错误描述

当我询问 AI 如何在 recipe 中引用环境变量时，得到的回答是使用 `${{ VAR_NAME }}`。这也是不正确的， 同样得通过 minijinja 模板引擎实现，语法如下：

```yaml
context:
  platform: {{ env.get("TARGET_PLATFORM", default="linux‑64") }}
```

### 4. pixi build 只能指定 platform，版本控制需另寻他法

我需要在一台机器上交叉处理其他平台的可执行文件，同时也需要指定处理的软件版本，再同时，我还希望使用 `pixi build` 来构建， 而不是 `rattler-build build`，这又踩了个坑， `pixi build` 目前是实验性质功能，因此它尚不能对 `rattler-build` 传递太多参数，目前只能通过环境变量来让 `pixi` 和 `rattler-build build` 获得一致的参数：

```bash
TARGET_PLATFORM=linux‑64 VERSION=1.1.35 pixi build
```

然后在 recipe 中通过 `{{ env.get("VERSION") }}` 来获取。

### 5. pixi 内置脚本的变量传递限制

由于要使用 `rattler-build bump-recipe`，为了让 action 和手动构建时输入的内容一致，也为了简单，我设置了 tasks，于是又是一坑，pixi 可以给 task 依赖传参数，但是不能指定环境变量（我前面要用环境变量来设置构建参数），当然这个问题还好，都是 shell，传不了我当参数传过去， 然后手动设置环境变量就是：

```toml
[tasks.bump]
cmd = "TARGET_PLATFORM={{ platform }} rattler-build bump-recipe"
args = [    
    { arg = "platform", default = "linux-64" },
]

[tasks.build]
cmd = "TARGET_PLATFORM={{ platform }} pixi build -t {{ platform }}"
args = ["platform"]
depends-on = [{ task = "bump", args = ["{{ platform }}"] }]
```

这里有要注意的一点是，build任务中必须环境变量和参数都指定平台信息，`pixi build` 似乎并不读`TARGET_PLATFORM`这个环境变量，因此要`-t`来指定平台，如果不指定，实际上会使用当前平台作为默认值，导致跨平台构建失败。

### 6. prefix.dev 缺少查询 API，需用 pixi search 解析

目前 prefix.dev 似乎没有提供专门的 API 来指定 conda 包查询对应版本信息。这个好解决，解析 `pixi search` 命令的输出就成，shell 作为万能胶水还是必不可少啊...

```bash
version=$(pixi search -q --no‑progress -p linux‑64 -c https://prefix.dev/sylens opencode | grep -oP '\K[0‑9]+\.[0‑9]+\.[0‑9]+' | head -1)
```

### 7. bump‑recipe 的重复下载问题

`rattler‑build bump‑recipe` 在更新 sha256 时会下载一次源码包，而随后的 `pixi build` 又会再次下载相同的包。这种重复下载在网速较慢或包较大时会显著增加构建时间。不过，这个目前确实无解，还好不影响功能实现。

## 完整实现

经过多次调整，最终的 pixi 项目配置如下：

```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda‑forge"]
name = "opencode‑recipe"
platforms = ["linux‑aarch64", "linux‑64"]
version = "0.1.0"
preview = ["pixi‑build"]

[dependencies]
rattler‑build = "*"

[tasks.bump]
cmd = "export TARGET_PLATFORM={{ platform }} && rattler‑build bump‑recipe"
args = [    
    { arg = "platform", default = "linux‑64" },
]

[tasks.build]
cmd = "export TARGET_PLATFORM={{ platform }} && pixi build"
args = ["platform"]
depends‑on = [{ task = "bump", args = ["{{ platform }}"] }]

[package.build.backend]
name = "pixi‑build‑rattler‑build"
version = "0.3.*"
```

recipe 文件的关键部分：

```yaml
context:
  name: opencode
  version: {{ env.get("VERSION", default="1.1.35") }}
  platform: {{ env.get("TARGET_PLATFORM", default="linux‑64") }}

package:
  name: {{ name }}
  version: {{ version }}

source:
  url: https://github.com/anomalyco/opencode/releases/download/v{{ version }}/opencode‑{{ "linux‑x64" if platform == "linux‑64" else "linux‑arm64" }}.tar.gz
  sha256: {{ sha256 }}
  file_name: opencode.tar.gz

build:
  number: 0
  string: h{{ hash }}_{{ platform }}
  script: |
    tar xvf opencode.tar.gz
    mkdir -p $PREFIX/bin
    mv opencode $PREFIX/bin/opencode
    chmod 755 $PREFIX/bin/opencode
```

GitHub Action 的核心步骤：

```yaml
- name: Run Pixi Task
  run: |
    cd opencode
    
    # 获取当前已上传版本
    current_version=$(pixi search -q --no‑progress -p ${{ matrix.platform }} -c https://prefix.dev/sylens opencode | grep -oP '\K[0‑9]+\.[0‑9]+\.[0‑9]+' | head -1)
    
    # 获取上游最新版本
    upstream_version="${{ steps.origin_release.outputs.release }}"
    
    if [ "$upstream_version" != "v$current_version" ]; then
      echo "检测到新版本: $upstream_version (当前: $current_version)"
      export VERSION=${upstream_version#v}
      pixi run build ${{ matrix.platform }}
      pixi upload prefix opencode*.conda -c sylens
    else
      echo "已是最新版本"
    fi
```