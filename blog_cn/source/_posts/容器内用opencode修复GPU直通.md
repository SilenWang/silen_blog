---
title: 容器内用 opencode + deepseek-v4-flash 修复 GPU 直通
categories: Coding
date: 2026-08-04 16:00:00
tags:
  - opencode
  - deepseek
  - GPU 直通
  - bwrap
  - Docker
---

之前已经成功把claude-science装进了容器，基本的功能使用没有太大问题，但是想做一点分子对接之类的事情时，由于容器内缺少显卡，无法进行，但是等我真按照之前的方式，拉取nvidia的官方镜像来构建容器后，诡异的情况发生了，CS中的所有语言kernel，pyhton，perl，R全部不能用，整个运行环境直接废了。

这真的超出我的理解范围，在手动请教Gemini无果后，我决定让 **opencode + deepseek-v4-flash** 自己处理——不是我在外面远程敲命令指挥，而是 AI agent 就在容器里，自己看日志、自己二分参数、自己反编译解析器、自己改 entrypoint。这体验还挺奇妙的，记录一下。

<!-- more -->

先根据AI的描述补充下背景：claude-science 会通过 bwrap 把每个任务隔离进沙箱，bwrap本身就是一种轻量的容器/沙箱技术，因此我的使用环境，相当于是在容器内嵌套使用容器；bwrap 沙箱需要直通宿主（容器）的 GPU，也就是要把 `/dev/nvidia*` 和 nvidia-smi 注入进去。如果不在容器内，这应该很顺，但在我这个使用条件下，出现了两个层面的问题。

## 问题一：bwrap 报错 + PID-namespace DEGRADED

第一个现象是，bwrap 起沙箱时报 `Can't mount proc on /newroot/proc`，同时日志里出现 "PID-namespace DEGRADED" 警告。

一开始按惯例怀疑是 CDI（Container Device Interface）把 privileged 搞失效了，于是各种排列组合：

- 加 `security_opt: seccomp=unconfined` / `apparmor=unconfined`
- 加 `cap_add: SYS_ADMIN`
- 逐步二分 bwrap 参数
- 对比两个容器的 mount 表
- 手工 umount `/proc/driver/nvidia/params` 后再跑探针

前几招都没用，直到对比 mount 表才发现不对劲：NVIDIA toolkit 竟然在 `/proc` 内部挂了子挂载。

**根因**：NVIDIA toolkit 的 `/proc/driver/nvidia/params` 是一个 tmpfs 子挂载，bwrap 要为嵌套的 pid namespace 重挂 proc 时，遇到 `/proc` 下面的子挂载就挂不上去了——这就是 `Can't mount proc on /newroot/proc` 和 DEGRADED 警告的来源。换句话说，这是 NVIDIA 的锅，不是权限配置的问题。

**修复**：在 entrypoint 里新增 `cleanup_proc_mounts()`，启动之前先把 `/proc/*` 下面的子挂载全部卸载掉。实测 nvidia-smi 和 CUDA 都不受影响，bwrap 也能正常挂上 proc。

## 问题二：沙箱里还是没有 GPU

proc 挂载解决之后，以为万事大吉，结果沙箱里一跑 nvidia-smi，退出码 9，`/dev/nvidia*` 一个都没有。设备、权限、配置看起来全对，就是不见 GPU。

又是一轮排查：

- 检查容器内设备节点是否存在
- 手工模拟带 GPU 绑定的 bwrap 沙箱
- 从二进制里反推配置 schema
- 设置 `gpu_enabled=true`（entrypoint 里走 `apply_gpu_config` + `ENABLE_GPU`）
- 把 `require_token=false` 想直接查 API
- 试 `NVIDIA_VISIBLE_DEVICES=void` 环境变量
- 最后只能反编译 sessionGpu 解析器（XHG）

折腾了一圈，终于在一段解析器逻辑里看到了真相：

```
sessionGpu = gpu_enabled && (gpu_mode != "off")
```

简单来说，即使容器配置、设备节点、权限全对，只要这个会话创建在没有 GPU 的时候，它的 gpu_mode 会保持 off，有了GPU后，claude-science 也不会把 GPU 注入沙箱。因此最简单的方式就是重新开会画。

当然 AI 还是找到了修复现有会画的方式：直接改数据库，把所有根会话的 `gpu_mode: off → on`。

## 用 opencode 直接在容器里调试的感受

这次的整个过程，是 opencode + deepseek-v4-flash 直接在容器内完成的：它自己跑命令复现、自己二分 bwrap 参数、自己对比 mount 表、自己反编译二进制、自己改 entrypoint 和 DB 迁移。我其实只负责检查到底能不能用，和把修改提交保存。

几点体会：

- 这类"环境类"问题特别适合 agent 在容器里就地调试，因为它能真实复现、实时验证，不用靠人来来回回贴日志。
- deepseek-v4-flash 在这种长链路排查里表现不错，尤其是"从现象反推配置/代码"的那几步。
- 当然，还得是在准备好的隔离环境中操作，这次它最终是顺利解决问题了，单保不齐以后那天会变成直接搞崩系统。。。
