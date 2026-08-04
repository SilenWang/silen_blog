---
title: Fixing GPU Passthrough in a Container with opencode + deepseek-v4-flash
categories: Coding
date: 2026-08-04 16:00:00
tags:
  - opencode
  - deepseek
  - GPU passthrough
  - bwrap
  - Docker
---

I recently worked on the sandbox environment for claude-science running in a container. The goal was simple: let tasks running inside the bwrap sandbox use the GPU. The whole process of investigation, diagnosis, and fixing was carried out entirely inside the container with **opencode + deepseek-v4-flash**—not me remotely typing commands and directing it, but an AI agent living inside the container that read the logs, bisected parameters, decompiled the parser, and edited the entrypoint itself. It was a strange and interesting experience; here's a record of it.

<!-- more -->

First, some background. claude-science isolates each task into a sandbox via bwrap, and the sandbox needs the host (container) GPU passed through, meaning `/dev/nvidia*` and nvidia-smi must be injected into it. In an ideal world this should be straightforward, but two completely independent layers of problems surfaced, and each one was buried quite deep.

## Problem 1: bwrap error + PID-namespace DEGRADED

The first symptom was bwrap failing to start the sandbox with `Can't mount proc on /newroot/proc`, along with a "PID-namespace DEGRADED" warning in the logs.

Following the usual instinct, I first suspected that CDI (Container Device Interface) had nullified privileged mode, so I tried every combination:

- Adding `security_opt: seccomp=unconfined` / `apparmor=unconfined`
- Adding `cap_add: SYS_ADMIN`
- Bisecting the bwrap parameters step by step
- Comparing the mount tables of two containers
- Manually unmounting `/proc/driver/nvidia/params` before running the probe again

The first few attempts were useless, until comparing the mount tables revealed something odd: the NVIDIA toolkit had mounted a sub-mount inside `/proc`.

**Root cause**: NVIDIA toolkit's `/proc/driver/nvidia/params` is a tmpfs sub-mount. When bwrap tried to remount proc for a nested PID namespace, the sub-mount under `/proc` blocked it—that's the source of both `Can't mount proc on /newroot/proc` and the DEGRADED warning. In other words, it was NVIDIA's fault, not a permission configuration problem.

**Fix**: Added `cleanup_proc_mounts()` to the entrypoint, which unmounts all sub-mounts under `/proc/*` before startup. In practice, nvidia-smi and CUDA are unaffected, and bwrap can now mount proc normally.

## Problem 2: Still no GPU inside the sandbox

After fixing the proc mount, I thought it was all done. But when I ran nvidia-smi inside the sandbox, it exited with code 9, and there were no `/dev/nvidia*` nodes at all. The devices, permissions, and configuration all looked correct—the GPU just wasn't there.

Another round of investigation:

- Checked whether the device nodes exist in the container
- Manually simulated a bwrap sandbox with GPU bindings
- Reverse-engineered the config schema from the binary
- Set `gpu_enabled=true` (via `apply_gpu_config` + `ENABLE_GPU` in the entrypoint)
- Set `require_token=false` to query the API directly
- Tried the `NVIDIA_VISIBLE_DEVICES=void` environment variable
- Finally resorted to decompiling the sessionGpu parser (XHG)

After a lot of digging, the truth finally appeared in a snippet of parser logic:

```
sessionGpu = gpu_enabled && (gpu_mode != "off")
```

**Root cause**: Our session was created at 16:58, at which point the GPU hadn't been enabled yet. When the session was created, `_original_input.gpu_mode="off"` was persisted into the root session's DB, and this session-level value takes precedence over the global `gpu_enabled=true`. The parser then worked its way through the logic and switched off GPU passthrough for the sandbox—the sandbox wasn't even given a `--dev-bind`.

In other words: even if the container configuration, device nodes, and permissions are all correct, as long as that session's gpu_mode is still off, claude-science will not inject the GPU into the sandbox. That's the real reason behind "config enabled but GPU still invisible".

**Fix**: Migrated the DB, flipping `gpu_mode: off → on` for all root sessions; meanwhile, `apply_gpu_config()` in the entrypoint now runs this migration automatically before startup, and new sessions default to on.

## Summary of the two root causes

Looking back, there were two completely independent blockers affecting GPU usage:

1. **NVIDIA toolkit's `/proc/driver/nvidia/params` tmpfs sub-mount** → broke bwrap's nested PID namespace, showing up as `Can't mount proc` + DEGRADED warning.
2. **The session-level gpu_mode switch** → the parser `sessionGpu = gpu_enabled && (gpu_mode != "off")`, the session was created before the GPU was enabled, `gpu_mode="off"` got persisted and took precedence over the global config, so the sandbox never injected `/dev/nvidia*`.

One problem lives at the system layer (mounts), the other at the business logic layer (config persistence). They're completely unrelated, yet both had to be fixed for the GPU to actually work. The worst thing about debugging problems like this is when "everything looks configured correctly"—had opencode not been able to dig through the mount tables and decompile the parser directly inside the container, I'd probably have spent days guessing one tiny piece at a time.

## What it's like debugging with opencode right inside the container

The entire process this time was done by opencode + deepseek-v4-flash directly inside the container: it reproduced the issue by running commands, bisected the bwrap parameters, compared mount tables, decompiled the binary, and edited the entrypoint and DB migration itself. My role was basically providing context and confirming the direction.

A few takeaways:

- "Environment-type" problems are a great fit for an agent debugging in place inside the container, because it can reproduce the issue for real and verify fixes live, without humans going back and forth pasting logs.
- deepseek-v4-flash performed well in this kind of long-chain investigation, especially on the "work backwards from symptom to config/code" steps.
- But an agent's debugging still needs method: bisect, compare, reverse-engineer. These methodologies don't stop working just because the tool changed. AI just speeds up the execution.
