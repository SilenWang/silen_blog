---
title: 使用 pixi 和 rattler_build 打包AI帮我写的程序
categories: Coding
date: 2026-01-25 20:50:24
tags: [pixi, rattler-build, conda-build, Go, devssh]
---

之前我已经尝试了修改和添加recipe到conda的频道，这次试试把我的[DevSSH](https://github.com/SilenWang/DevSSH)打包上传到我自己的频道内。这次我想试试自己打conda包。

<!-- more -->

## 为什么选择pixi和rattler-build？

这次其实是个巧合... 本来想注册anaconda帐号用`conda-build`，结果他们官网不知道抽什么风，就是注册不成功，于是直接转prefix.dev了，正好学习一下更新的共苦

## rattler-build简介
`rattler-build`是开发`pixi`团队开发的新一代conda包构建工具，它用Rust编写，自然就拥有更快的构建速度，它一定程度上与conda-forge的构建系统兼容，并且与`pixi`互相进行了一定程度的集成。

## 配置pixi.toml

`pixi`的基本功能是配置开发环境，但是其实它也可以用来配置软件包的构建环境。进行相关的配置需要用到package关键字，同时需要在`[workspace]`内开启实验性质的`pixi-build`特性。

```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda-forge"]
name = "DevSSH"
platforms = ["linux-aarch64", "linux-64"]
version = "0.1.1"
preview = ["pixi-build"] # 需要开启pixi-build

[tasks]
build = {cmd = "go build -o bin/devssh cmd/devssh/main.go", cwd = "./"}

[activation.env]
CGO_ENABLED = "0"

[dependencies]
go = ">=1.25.4,<2"

# 下面是构建包的配置部分
[package]
name = "DevSSH"
version = "0.1.1"

[package.build.backend]
name = "pixi-build-rattler-build" # 构建的是go程序，没有专门构建器，所以使用rattler-build
version = "0.3.*" 
```

## 编写rattler-build的recipe

rattler-build使用YAML格式的recipe文件来定义如何构建包，它与`conda-build`下的recipe是兼容的，不过并不支持所有的`conda-build` recipe语法。

以下是为DevSSH编写的recipe：

```yaml
package:
  name: devssh
  version: 0.1.1

source:
  path: .
  use_gitignore: true 

build:
  number: 0
  script: |
    # 构建二进制文件
    go build -o "$PREFIX/bin/devssh" cmd/devssh/main.go

requirements:
  build:
    - ${{ compiler('go-nocgo') }} # 我的程序纯go没有c依赖，因此选择nocgo编译器
    - patchelf  # 执行post-process必须

about:
  homepage: https://github.com/SilenWang/DevSSH
  license: MPL-2.0
  license_file: LICENSE
  summary: 'A CLI tool to quickly set up remote development tools over SSH'
  description: |
    A CLI tool to quickly set up remote development tools over SSH
  repository: https://github.com/SilenWang/DevSSH
```

## 构建和上传包

### 1. 本地构建

如果配置文件设置无问题，执行构建命令应该就能得到conda格式的包了。

```bash
# 使用pixi构建
pixi build

# 指定目标平台构建
pixi build -t linux-64
```

如果有问题，就一点点问AI，不过值得注意的是，现在rattler-build可能不是太主流，得到的答案经常会是按照conda-build来答的，还是需要自己对着文档去鉴别尝试。

### 2. 配置prefix.dev账户并上传

上传前还是需要登录了，需要在账户下创建频道和设置API_TOKEN。

```bash
# 设置API令牌
pixi auth login --token YOUR_API_TOKEN prefix.dev
```

设置好了，直接把编译好的内容上传就是，网速好的话很快就完成了。

```bash
pixi upload prefix devssh-0.1.1-hb0f4dca_0.conda -c sylens
```

### 3. 通过频道安装包

用 pixi 或者 conda 都可以安装。

```bash
pixi global install -c https://prefix.dev/sylens devssh

conda install -c https://prefix.dev/sylens devssh
```
