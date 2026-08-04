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

I had already managed to get claude-science running in a container, and basic functionality worked without major issues. But when I wanted to do things like molecular docking, I couldn't—the container simply had no GPU. However, when I pulled NVIDIA's official image to build the container the same way I had before, something bizarre happened: every language kernel in CS—Python, Perl, R—stopped working, and the entire runtime environment was rendered useless.

This was truly beyond my understanding. After consulting Gemini manually with no luck, I decided to let **opencode + deepseek-v4-flash** handle it on its own—not me remotely typing commands and directing it, but an AI agent living inside the container that read the logs, bisected parameters, decompiled the parser, and edited the entrypoint itself. It was a strange and interesting experience; here's a record of it.

<!-- more -->

Let me first add some background based on the AI's description: claude-science isolates each task into a sandbox via bwrap, and bwrap itself is a lightweight container/sandbox technology. So my setup is essentially nesting containers inside a container. The bwrap sandbox needs the host (container) GPU passed through, meaning `/dev/nvidia*` and nvidia-smi must be injected into it. If it weren't inside a container, this should be straightforward, but under my usage conditions, problems appeared on two levels.

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

Put simply, even if the container configuration, device nodes, and permissions are all correct, as long as the session was created at a time when there was no GPU, its gpu_mode stays off—and once a GPU becomes available, claude-science still won't inject the GPU into the sandbox. So the simplest way out is to just open a new session.

Of course, the AI did find a way to fix the existing sessions too: modifying the database directly, flipping `gpu_mode: off → on` for all root sessions.

## What it's like debugging with opencode right inside the container

The entire process this time was done by opencode + deepseek-v4-flash directly inside the container: it reproduced the issue by running commands, bisected the bwrap parameters, compared mount tables, decompiled the binary, and edited the entrypoint and DB migration itself. All I actually did was check whether it worked and commit/save the changes.

A few takeaways:

- "Environment-type" problems are a great fit for an agent debugging in place inside the container, because it can reproduce the issue for real and verify fixes live, without humans going back and forth pasting logs.
- deepseek-v4-flash performed well in this kind of long-chain investigation, especially on the "work backwards from symptom to config/code" steps.
- Of course, it also has to operate in a prepared, isolated environment. This time it eventually solved the problem cleanly, but who knows—one day it might end up directly breaking the system...
