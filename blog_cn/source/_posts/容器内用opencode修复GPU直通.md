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

最近在容器里折腾 claude-science 的沙箱环境，目标很明确：让跑在 bwrap 沙箱里的任务能用上 GPU。整个排查、定位、修复的过程，全部是在容器内部，用 **opencode + deepseek-v4-flash** 直接完成的——不是我在外面远程敲命令指挥，而是 AI agent 就在容器里，自己看日志、自己二分参数、自己反编译解析器、自己改 entrypoint。这体验还挺奇妙的，记录一下。

<!-- more -->

先说下背景。claude-science 会通过 bwrap 把每个任务隔离进沙箱，沙箱需要直通宿主（容器）的 GPU，也就是要把 `/dev/nvidia*` 和 nvidia-smi 注入进去。理想状态下这应该很顺，但实际上冒出来两层完全独立的问题，而且每一个都藏得挺深。

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

**根因**：我们的会话创建于 16:58，而当时 GPU 还没启用。创建会话时，`_original_input.gpu_mode="off"` 被持久化到了根会话的 DB 里，并且这个会话级的值优先于全局的 `gpu_enabled=true`。于是解析器一路算下来，把沙箱的 GPU 直通关掉了——沙箱连 `--dev-bind` 都没加。

换句话说：即使容器配置、设备节点、权限全对，只要这个会话的 gpu_mode 还是 off，claude-science 就不会把 GPU 注入沙箱。这正是"配置已启用却仍看不到 GPU"的真正原因。

**修复**：对 DB 做迁移，把所有根会话的 `gpu_mode: off → on`；同时 entrypoint 的 `apply_gpu_config()` 现在会在启动前自动执行这个迁移，新会话默认就是 on。

## 两个根因小结

回头总结，真正影响 GPU 使用的是两个完全独立的阻塞点：

1. **NVIDIA toolkit 的 `/proc/driver/nvidia/params` tmpfs 子挂载** → 破坏了 bwrap 的嵌套 pidns，表现为 `Can't mount proc` + DEGRADED 警告。
2. **会话级 gpu_mode 开关** → 解析器 `sessionGpu = gpu_enabled && (gpu_mode != "off")`，会话创建于 GPU 启用之前，`gpu_mode="off"` 被持久化并优先于全局配置，导致沙箱完全不注入 `/dev/nvidia*`。

两个问题一个在系统层面（mount），一个在业务逻辑层面（配置持久化），完全不相干，却都必须修掉才能让 GPU 真正可用。调试这类问题最怕的就是"配置看起来都对"，这次要不是 opencode 能直接在容器里翻 mount 表、反编译解析器，光靠人肉一点一点猜，估计得耗上好几天。

## 用 opencode 直接在容器里调试的感受

这次的整个过程，是 opencode + deepseek-v4-flash 直接在容器内完成的：它自己跑命令复现、自己二分 bwrap 参数、自己对比 mount 表、自己反编译二进制、自己改 entrypoint 和 DB 迁移。我的角色基本就是提供上下文和确认方向。

几点体会：

- 这类"环境类"问题特别适合 agent 在容器里就地调试，因为它能真实复现、实时验证，不用靠人来来回回贴日志。
- deepseek-v4-flash 在这种长链路排查里表现不错，尤其是"从现象反推配置/代码"的那几步。
- 但 agent 排查也得有章法：二分、对比、反推，这些方法论不会因为换了工具就失效。AI 只是把执行提速了。
