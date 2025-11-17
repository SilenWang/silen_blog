---
title: Half Success - Compiling and Running Podman on FydeOS Using Pixi
categories: Others
date: 2025-11-11 01:30:59
tags: ['pixi', 'podman', 'FydeOS']
---

Recently, I attempted to use pixi to compile and run podman on FydeOS. While I made some progress, it ultimately only achieved half success. In this blog post, I'll share the entire process and the issues encountered.

<!-- more -->

## Background and Motivation

I use Fydetab Duo for local development and wanted to improve development efficiency through devpod. However, the Linux environment that comes with FydeOS has stability issues: after the device wakes from sleep, it often cannot access the Linux environment via the network, which severely impacts the development workflow.

As an alternative, I considered running podman outside the virtual machine to manage the development environment. This way, even if the Linux environment becomes unavailable, the development environment could still be maintained through podman. However, the podman package provided by conda-forge doesn't have an ARM64 version, making it impossible to use directly on Fydetab Duo (which is based on ARM architecture).

Therefore, I decided to clone the podman package build code from conda-forge and use the pixi environment management tool to compile the ARM64 version of podman myself. Pixi is a cross-platform package and environment manager that can help manage complex compilation dependencies.

## Environment Configuration and Compilation Process

First, I cloned the podman package build code from conda-forge and created a pixi project to manage the compilation environment. The `pixi.toml` file was configured as follows:

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

This configuration specifies project dependencies and system requirements, particularly for ARM64 architecture (linux-aarch64) and glibc version 2.39. The version specification in the `[system-requirements]` section is actually Pixi's virtual package specification mechanism, used to inform the dependency calculator about the current system's library versions during dependency calculation (everything depends on libc). The `channels` section specifies a local path because `podman` has another necessary dependency `netavark` that also doesn't have an Arm64 version, which needs to be compiled first to form a local Conda package before compiling `podman`.

## Compilation Challenges and Solutions

The conda package build code for the two software packages are located at `https://github.com/conda-forge/netavark-feedstock` and `https://github.com/conda-forge/podman-feedstock` respectively. After cloning them, compile them in sequence. When compiling both `netavark` and `podman`, there were dependency issues. We need to modify their `recipe/conda_build_config.yaml` files. The specific differences are shown below, essentially specifying to use sysroot as the C standard library.

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

After successful compilation, you can use `pixi add` to add podman to the environment, and then you'll have the podman command line available.

## Core Issue Leading to Runtime Failure

Everything went quite smoothly up to this point. However, when running `podman pull alpine`, it would prompt that `newgidmap` and `newuidmap` were missing. These are used for user and group ID mapping. This didn't stump me, as I'm used to finding workarounds - I just copied them from the Linux virtual machine. However, even with these programs present, I still encountered the following error:
```
newuidmap: Could not set caps
cannot set up namespace using "/usr/bin/newuidmap": should have setuid or have filecaps setuid
```

At this point, I encountered an insurmountable obstacle: **FydeOS/ChromeOS is designed from the ground up to not allow UID/GID mapping in user namespaces**.

I tried setting the setuid bit for `newuidmap` and `newgidmap` using `sudo chmod u+s` as suggested by AI, but the system still refused to execute them. The kernel log showed:
```
SafeSetID: Operation requires CAP_SETUID, which is not available to UID 1000 for operations besides approved set*uid transitions
```

Here's what AI told me about this:

> 1. SafeSetID is a ChromeOS-specific security mechanism (LSM/kernel patch) used to restrict setuid / UID mapping
> 2. CAP_SETUID is unavailable → Even if newuidmap has setuid or setcap, regular users (UID 1000) cannot modify UID mapping in rootless Podman
> 3. Other logs are normal remount or USB device information, unrelated to rootless Podman

Since I don't have enough knowledge to verify the accuracy of this information, I could only test by compiling Podman in the same way within the Linux virtual machine. The test confirmed that Podman in the virtual machine runs normally, which at least proves that there are indeed special restrictions on the host system.

## Conclusion

This attempt proved that using pixi to manage podman dependencies on FydeOS is feasible, and the compilation process can be completed successfully. However, due to deep kernel security restrictions in FydeOS (particularly the SafeSetID module), the runtime environment cannot meet podman's requirements for user namespaces, so this can only be considered "half success."

Therefore, to achieve my goal, I'll need to compile the Openfyde image myself, directly modifying parts of the kernel during compilation and including podman directly.