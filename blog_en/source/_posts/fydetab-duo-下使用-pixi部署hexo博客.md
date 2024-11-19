---
title: Deploying Hexo Blog with Pixi on Fydetab Duo
categories: Script
date: 2024-08-04 13:04:05
tags: ['hexo', 'pixi', 'fydeos', 'fydetab duo']
---


I've had the Fydetab Duo for over a month, and I'm still trying to find ways to make the most use of this device. At the very least, I'd like to write some blogs during my downtime with it. Hence the following exploration.

<!-- more -->

To be honest, as the Chinese version of Chrome OS, Fyde OS, the recommended method for using Linux should be directly through the Linux subsystem. However, within this month, the Linux subsystem seems to still have some issues that affect usage, such as not being able to connect to the internet for no apparent reason after waking up, mouse display anomalies within the subsystem, and relatively slow subsystem startup, which are minor but impactful issues. Plus, compared to Chrome OS, the Fyde OS team still has more lenient restrictions on the host system (at least sudo does not require switching to a different tty), so I tried to deploy the tools I needed on the host system.

## Installing Pixi

- To install `pixi`, additional environment variable configuration is required because under the Chromium OS-based system, not all places can store and run executables. `/usr/local/share` is the only choice, so you need to add the following content to the environment variable file, and then run `curl -fsSL https://pixi.sh/install.sh | bash`

```bash
export PATH=/usr/local/share/pixi/bin:$PATH
export PIXI_HOME=/usr/local/share/pixi
```

## Installing code-server and git

The original choice was Microsoft's official `vscode-cli` tool, but possibly due to the differences in file structure and environment configuration between Fyde OS and general Linux distributions, it keeps getting stuck at the step of downloading the vscode server when connecting to the service, making it unusable in practice. Therefore, I temporarily used pixi to install `code-server`.

## Special Configuration for ARM Platforms

- Fydeos Tab is an ARM architecture device, so the configuration for `pixi` needs to include the `linux-aarch64` platform. In practice, installing yarn has resulted in unrecognized sub-environment python issues, so Fydetab Duo does not use yarn and reverts to npm. Also, since Fyde OS does not have `make`, `gcc`, `gxx` by default, these components need to be added to the `linux-aarch64` dependencies.
- The system information files in Fyde OS have minor differences from conventional Linux, which can cause hexo's module to incorrectly parse system information, rendering `hexo -v` and the `hexo g` that calls it inoperative. Therefore, it is necessary to add a fixing statement in the deployment task.

```toml
[target.linux-aarch64.tasks]
build_blog = {cmd = "npm install -g hexo-cli && npm install --python=$PWD/.pixi/envs/default/bin/python2 && sed -i '34s/distro\\[1\\]/distro/' $PWD/.pixi/envs/default/lib/node_modules/hexo-cli/dist/console/version.js", cwd = "."}
```

## Cloning And Deploying the Blog

These are the additional configurations needed. Simply clone the project and deploy it directly, and you can continue writing blogs.

```bash
git clone https://github.com/SilenWang/silen_blog
pixi run build_blog
```
