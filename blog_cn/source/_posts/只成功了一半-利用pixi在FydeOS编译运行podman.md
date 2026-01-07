---
title: 只成功了一半----利用pixi在FydeOS编译运行podman
categories: Others
date: 2025-11-11 01:30:59
tags: ['pixi', 'podman', 'FydeOS']
---

最近在FydeOS上尝试使用pixi来编译和运行podman，虽然取得了一些进展，但最终只成功了一半。在这篇博客中，我将分享整个过程以及遇到的问题。

<!-- more -->

## 背景与动机

我使用Fydetab Duo进行本地开发，希望通过devpod来提高开发效率。然而，FydeOS自带的Linux环境存在稳定性问题：设备从休眠状态唤醒后，经常无法通过网络访问Linux环境，这严重影响了开发工作流程。

作为替代方案，我考虑在虚拟机外运行podman来管理开发环境。这样即使Linux环境不可用，也能通过podman维持开发环境的运行。但是，conda-forge提供的podman包没有ARM64架构版本，无法直接在Fydetab Duo（基于ARM架构）上使用。

因此，我决定克隆conda-forge的podman包构建代码，使用pixi环境管理工具自行编译ARM64版本的podman。Pixi是一个跨平台的包管理器和环境管理器，能够帮助管理复杂的编译依赖关系。

## 环境配置与编译过程

首先，我克隆了conda-forge的podman包构建代码，并创建了一个pixi项目来管理编译环境。配置`pixi.toml`文件如下：

```toml
[workspace]
authors = ["Sylens <qiumin14@163.com>"]
channels = ["conda-forge", "file:///usr/local/share/project/podman/.pixi/envs/default/conda-bld"]
name = "podman"
platforms = ["linux-aarch64"]
version = "0.1.0"

[tasks]

[system-requirements]
libc   = { family = "glibc", version = "2.39" }

[dependencies]
boa = ">=0.17.0,<0.18"
mamba = ">=1.5.12,<2"
conda-build = ">=24.5.1,<25"
diffutils = ">=3.12,<4"
podman = ">=5.6.2,<6"
libcap = ">=2.76,<3"
```

这个配置指定了项目依赖和系统要求，特别是针对ARM64架构（linux-aarch64）和glibc 2.39版本。`[system-requirements]`部分的版本指定，其实是Pixi的一种虚拟包指定机制，用于在计算依赖时，告诉依赖计算器当前系统的系统库版本（万物依赖libc）。而`channels`中指定了一个本地的路径，因为`podman`还有一个必要依赖`netavark`也没有Arm64版本，需要先行编译形成本地Conda包后，才能开始编译`podman`。

## 编译解决方案

两个软件的conda包构建代码分别在`https://github.com/conda-forge/netavark-feedstock`和`https://github.com/conda-forge/podman-feedstock`，克隆下来后依次编译。在编译`netavark`和`podman`时，都会有依赖问题，我们需要修改他们的`recipe/conda_build_config.yaml`文件，具体差异别见下，其实就是指定使用sysroot作为C标准库。

```
index d402d1d..7dbfd74 100644
--- a/recipe/conda_build_config.yaml
+++ b/recipe/conda_build_config.yaml
@@ -1,2 +1,5 @@
+c_stdlib:
+  - sysroot
+
 c_stdlib_version:  # [linux64]
   - 2.17           # [linux64]
```

编译成功后，就可以使用`pixi add`添加podman到环境，就有podman命令行了

## 运行失败的核心问题

事情到这都还挺顺利的，然而如果运行`podman pull alpine`，会提示没有`newgidmap`程序和`newuidmap`。他们是用来做用户和组编号映射的。这难不倒常年借东墙补西墙的我，Linux虚拟机内复制一个出来就行。然而，在有这两个程序后，依然会看到下面的错误：
```
newuidmap: Could not set caps
cannot set up namespace using "/usr/bin/newuidmap": should have setuid or have filecaps setuid
```

到此位置，遇到了本次无法克服的障碍，**FydeOS/Chromeos的从设计上，就不允许用户命名空间中的UID/GID映射**。

我尝试根据AI提示，通过`sudo chmod u+s`为`newuidmap`和`newgidmap`设置了setuid位，系统仍然拒绝执行。内核日志显示：
```
SafeSetID: Operation requires CAP_SETUID, which is not available to UID 1000 for operations besides approved set*uid transitions
```

以下是AI给我的回答：

>1.SafeSetID 是 ChromeOS 特有的安全机制（LSM/内核补丁），用于 限制 setuid / UID 映射
>2.CAP_SETUID 不可用 → 即使 newuidmap 有 setuid 或 setcap，普通用户（UID 1000）也无法在 rootless Podman 中修改 UID 映射
>3.其他日志都是正常 remount 或 USB 设备信息，与 rootless Podman 无关

由于我没有足够的知识来判断它说的真假，所以我只能进Linux虚拟机以一样的方式编译Podman试了试，实测，虚拟机内的Podman是能正常运行的，这至少说明，宿主机部分确实有特别的限制。

## 总结

这次尝试证明使用pixi在FydeOS上管理podman依赖是可行的，编译过程也能成功完成。但由于FydeOS深层的内核安全限制（特别是SafeSetID模块），运行时环境无法满足podman对用户命名空间的要求，因此只能算"成功了一半"。

所以，如果要实现我的目标，我只能自行去编译Openfyde镜像，直接在编译的时候修改内核的部分，并且直接加入podman了。
