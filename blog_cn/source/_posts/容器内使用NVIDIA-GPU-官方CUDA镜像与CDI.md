---
title: 容器内使用 NVIDIA GPU 的两种方式：官方 CUDA 镜像 vs CDI
categories: Coding
date: 2026-08-05 22:00:00
tags:
  - NVIDIA
  - CUDA
  - CDI
  - Docker
  - GPU
  - 容器技术
---

设置 Cluade Science 容器的时候，了解到了新的容器内使用 N 社 GPU的方式，记录一下（当然大部分都是AI写的，作为我自己检索查看的材料用）。

<!-- more -->

在容器里用上 NVIDIA GPU，早年间最常见的做法，就是直接拉 NVIDIA 官方、自带 CUDA 的容器镜像（`nvidia/cuda` 系列），配合 NVIDIA Container Toolkit 来跑。而这两年 CDI（Container Device Interface）逐渐成了更"现代"的方式，Docker、以及多种容器技术都在积极的原生支持它。

两种方式表面上看都是"容器里能用 GPU"，但底层的技术实现和日常使用的体验，差别其实非常大。

## 老方式：官方 CUDA 镜像 + NVIDIA Container Runtime

先说老方式。NVIDIA 官方在 Docker Hub 上维护着一堆 `nvidia/cuda` 镜像，比如 `nvidia/cuda:12.4.1-devel-ubuntu22.04`。这些镜像把 CUDA 工具链、运行时库（`libcudart`、`libcublas` 这些）直接打包进了镜像里，拉下来就能编译、运行 CUDA 程序。

但光有镜像还不够——镜像里没有的驱动部分（`libcuda.so`、设备节点 `/dev/nvidia*`）需要从宿主机注入。这就要装 NVIDIA Container Toolkit（早年叫 `nvidia-docker2`），并把它配置成 Docker 的运行时：

```bash
# 安装 toolkit 后配置 docker 运行时
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# 跑一个 CUDA 容器
docker run --rm --gpus all \
  -e NVIDIA_VISIBLE_DEVICES=all \
  nvidia/cuda:12.4.1-base-ubuntu22.04 \
  nvidia-smi
```

### 实现原理

老方式的核心是 **hook 机制**。整个链路大致是这样的：

1. `nvidia-container-runtime` 作为 Docker 的运行时被调用，它本质是 `runC` 的薄包装；
2. 运行时往 OCI 规范（`config.json`）里注入一个 `prestart` hook（对应 `/usr/share/containers/oci/hooks.d/oci-nvidia-hook.json`）；
3. hook 调用 `nvidia-container-cli`，根据 `NVIDIA_VISIBLE_DEVICES` 环境变量决定注入哪些设备；
4. CLI 把宿主机上对应的设备节点、驱动库（`libcuda.so` 等）、以及 `nvidia-smi` 等工具挂载/复制进容器。

所以这个模式下有几个关键点：

- **设备选择靠环境变量**：`NVIDIA_VISIBLE_DEVICES` 可以写 `all`、`0`、`1`、UUID 等；
- **CUDA 在镜像里，驱动在宿主机**：镜像打包的是 CUDA 用户态库，驱动由宿主机通过 hook 注入，两者在运行时"拼"起来；
- **跟运行时强绑定**：这套 hook 机制是 NVIDIA 私有的，Docker 要配置 `nvidia` 运行时，Podman 默认也能识别 hook，但每个运行时都得单独适配。

### 老方式的痛点

用久了会发现一些问题：

- **镜像巨大**：CUDA 工具链打进去，动辄几个 GB，拉取和存储都是负担；
- **版本耦合**：镜像里的 CUDA 版本和宿主机驱动版本需要匹配，升级驱动后老镜像可能就不能用了；
- **跟运行时绑定太紧**：hook 机制是 NVIDIA 一家之言，换别的运行时或厂商的设备就得重来；
- **rootless 支持差**：hook 机制在 rootless 容器场景下兼容性不好。

## 新方式：CDI（Container Device Interface）

CDI 是 CNCF 下的一个开放规范（`cncf-tags/container-device-interface`），想解决的就是上面那堆问题。它的思路是把"设备如何暴露给容器"这件事标准化，让运行时、容器引擎、编排器都能用同一套语言来描述设备。

### 基本概念

- 设备用**全限定名**标识，格式是 `vendor.com/class=unique_name`，NVIDIA 就是 `nvidia.com/gpu=0`、`nvidia.com/gpu=1:0`（MIG 设备）、`nvidia.com/gpu=all` 这种；
- 设备信息以规范文件（`.yaml`/`.json`）形式放在默认目录 `/etc/cdi`（静态配置）和 `/var/run/cdi`（动态更新）；
- **CDI 只负责让容器"感知"设备**，资源管理（调度、配额）是编排器的活，CDI 不管。

一个 CDI 规范文件大概长这样（节选）：

```yaml
cdiVersion: 0.6.0
kind: nvidia.com/gpu
devices:
  - name: "0"
    containerEdits:
      deviceNodes:
        - hostPath: /dev/nvidia0
          path: /dev/nvidia0
      mounts:
        - hostPath: /usr/lib/x86_64-linux-gnu/libcuda.so.1
          path: /usr/lib/x86_64-linux-gnu/libcuda.so.1
```

可以看到，CDI 把"注入哪些设备节点、挂载哪些库、设哪些环境变量"都声明成了数据。容器引擎读到设备名，查一下规范文件，就能自己把对应的内容拼进 OCI 规范，不需要再走 NVIDIA 私有的 hook。

### NVIDIA 如何生成 CDI 规范

NVIDIA Container Toolkit 从 **v1.12.0** 起支持生成 CDI 规范。从 **v1.18.0** 起，systemd 服务 `nvidia-cdi-refresh` 会在 Toolkit 安装/升级、GPU 驱动安装/升级、系统重启时，自动生成并更新位于 `/var/run/cdi/nvidia.yaml` 的规范：

```bash
# 查看当前可用的 CDI 设备
nvidia-ctk cdi list
```

有个注意点：`nvidia-cdi-refresh` 目前**不自动处理驱动移除和 MIG 设备重配置**这两种场景，需要手动重新生成：

```bash
sudo nvidia-ctk cdi generate --output=/var/run/cdi/nvidia.yaml
```

### 怎么用 CDI

**Docker** 从 **25.0.0** 起支持 CDI，**28.2.0** 起默认启用。中间的版本（25.0.0~28.1.1）需要在 `/etc/docker/daemon.json` 里手动开启：

```json
{
  "features": {
    "cdi": true
  }
}
```

运行时候直接：

```bash
docker run --rm \
  --device nvidia.com/gpu=all \
  ubuntu:24.04 nvidia-smi
```

注意这里跑的是普通的 `ubuntu` 镜像，CUDA 库由 CDI 规范从宿主机注入——这正好戳中老方式"镜像巨大"的痛点。

### 不原生支持 CDI 的运行时怎么办

对不原生支持 CDI 的运行时（比如没开 CDI 特性的 Docker），可以把 NVIDIA Container Runtime 配置成 `cdi` 模式，此时设备选择仍然用 `NVIDIA_VISIBLE_DEVICES` 环境变量，但注入走 CDI：

```bash
sudo nvidia-ctk config --in-place --set nvidia-container-runtime.mode=cdi

docker run --rm -ti --runtime=nvidia \
  -e NVIDIA_VISIBLE_DEVICES=nvidia.com/gpu=all \
  ubuntu nvidia-smi -L
```

一个需要记住的限制：**CDI 注入和 NVIDIA Container Runtime 的 hook 机制是互斥的**。如果存在 `/usr/share/containers/oci/hooks.d/oci-nvidia-hook.json`，就别同时用 CDI 注入，也不要再设 `NVIDIA_VISIBLE_DEVICES`，否则会打架。

## 两种方式的对比

| 维度 | 老方式：官方 CUDA 镜像 | 新方式：CDI |
| --- | --- | --- |
| 设备注入机制 | NVIDIA 私有 prestart hook + `nvidia-container-cli` | 开放规范的 CDI 文件，运行时直接解析 |
| 设备选择方式 | `NVIDIA_VISIBLE_DEVICES` 环境变量 | `--device nvidia.com/gpu=xxx` 或注解/CRI 字段 |
| CUDA 来源 | 打包进镜像，镜像巨大 | 宿主驱动 + CDI 注入，可用普通镜像 |
| 版本耦合 | 镜像 CUDA 与驱动需匹配，升级敏感 | 规范由 `nvidia-cdi-refresh` 自动跟随驱动 |
| 运行时兼容 | 依赖 NVIDIA hook，逐运行时适配 | Podman/Docker/containerd/CRI-O 原生统一支持 |
| rootless | 支持差 | 支持好 |
| 生态开放度 | NVIDIA 私有 | CNCF 开放规范，其他厂商设备也可用 |
| 资源管理 | 无标准 | 明确交给编排器（如 K8s device plugin） |

## 小结

老方式简单直接，很多人（包括我）第一套 GPU 容器就是这么跑起来的，现在依然能用，适合快速验证、个人使用。但如果你的环境要上 K8s、要支持 rootless、或者有多运行时并存，CDI 是明显更现代、更省心的方向。

不过话说回来，CDI 也不是银弹——`nvidia-cdi-refresh` 不覆盖驱动移除和 MIG 重配置，规范偶尔得手动刷；而且它要求运行时侧支持，Docker 老版本还得手动开特性开关。但总体趋势很明确：设备描述正在从"各家私有的 hook"走向"开放、声明式的规范"。
