---
title: Packaging an AI-Assisted Program with pixi and rattler-build
categories: Coding
date: 2026-01-25 20:50:24
tags: [pixi, rattler-build, conda-build, Go, devssh]
---

Previously, I had tried modifying and adding recipes to conda channels. This time, I want to package my [DevSSH](https://github.com/SilenWang/DevSSH) and upload it to my own channel. I decided to try building a conda package myself.

<!-- more -->

## Why choose pixi and rattler-build?

This was actually a coincidence... I originally wanted to register an Anaconda account to use `conda-build`, but their website was acting up and I couldn't complete the registration. So I switched directly to prefix.dev, which gave me a chance to learn about newer tools.

## Introduction to rattler-build
`rattler-build` is a next-generation conda package builder developed by the team behind `pixi`. Written in Rust, it naturally offers faster build speeds. It is somewhat compatible with conda-forge's build system and has a certain level of integration with `pixi`.

## Configuring pixi.toml

`pixi`'s core functionality is to set up development environments, but it can also be used to configure build environments for software packages. To do this, you need to use the `package` keyword and enable the experimental `pixi-build` feature within the `[workspace]` section.

```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda-forge"]
name = "DevSSH"
platforms = ["linux-aarch64", "linux-64"]
version = "0.1.1"
preview = ["pixi-build"] # Need to enable pixi-build

[tasks]
build = {cmd = "go build -o bin/devssh cmd/devssh/main.go", cwd = "./"}

[activation.env]
CGO_ENABLED = "0"

[dependencies]
go = ">=1.25.4,<2"

# Below is the configuration for building the package
[package]
name = "DevSSH"
version = "0.1.1"

[package.build.backend]
name = "pixi-build-rattler-build" # Building a Go program, no specialized builder, so use rattler-build
version = "0.3.*" 
```

## Writing a rattler-build recipe

rattler-build uses YAML‑format recipe files to define how a package is built. It is compatible with recipes from `conda-build`, though it does not support all `conda-build` recipe syntax.

Here is the recipe I wrote for DevSSH:

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
    # Build the binary
    go build -o "$PREFIX/bin/devssh" cmd/devssh/main.go

requirements:
  build:
    - ${{ compiler('go-nocgo') }} # My program is pure Go with no C dependencies, so choose the nocgo compiler
    - patchelf  # Required for post‑processing

about:
  homepage: https://github.com/SilenWang/DevSSH
  license: MPL-2.0
  license_file: LICENSE
  summary: 'A CLI tool to quickly set up remote development tools over SSH'
  description: |
    A CLI tool to quickly set up remote development tools over SSH
  repository: https://github.com/SilenWang/DevSSH
```

## Building and uploading the package

### 1. Local build

If the configuration files are set up correctly, running the build command should produce a conda‑format package.

```bash
# Build with pixi
pixi build

# Build for a specific target platform
pixi build -t linux-64
```

If you encounter issues, you can ask AI step by step. However, note that rattler-build is not yet mainstream, so answers you get may often be based on conda-build. You'll need to cross‑check with the documentation and experiment.

### 2. Setting up a prefix.dev account and uploading

Before uploading, you need to log in, create a channel under your account, and set up an API_TOKEN.

```bash
# Set the API token
pixi auth login --token YOUR_API_TOKEN prefix.dev
```

Once configured, you can directly upload the built package. With a good internet connection, it completes quickly.

```bash
pixi upload prefix devssh-0.1.1-hb0f4dca_0.conda -c sylens
```

### 3. Installing the package via the channel

You can install it with either pixi or conda.

```bash
pixi global install -c https://prefix.dev/sylens devssh

conda install -c https://prefix.dev/sylens devssh
```