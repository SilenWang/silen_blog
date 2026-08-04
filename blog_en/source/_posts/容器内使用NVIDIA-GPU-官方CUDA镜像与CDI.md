---
title: Two Ways to Use NVIDIA GPU in Containers: Official CUDA Images vs CDI
categories: Coding
date: 2026-08-05 22:00:00
tags:
  - NVIDIA
  - CUDA
  - CDI
  - Docker
  - Podman
  - GPU
  - Containers
---

Using an NVIDIA GPU inside containers has changed quite a bit over the past couple of years. In the early days, the most common approach was to pull NVIDIA's official CUDA container images (the `nvidia/cuda` series) and run them with the NVIDIA Container Toolkit. In the last two years, CDI (Container Device Interface) has become the more "modern" way, with native support in Podman, Docker, containerd, and CRI-O.

On the surface, both approaches let you "use a GPU in a container," but the underlying technical implementation and day-to-day experience are actually very different. This post is about my understanding, focusing on the differences in **implementation** and in **features and usability** between the two.

<!-- more -->

## The old way: official CUDA images + NVIDIA Container Runtime

Let's start with the old way. NVIDIA maintains a bunch of `nvidia/cuda` images on Docker Hub, such as `nvidia/cuda:12.4.1-devel-ubuntu22.04`. These images bake the CUDA toolchain and runtime libraries (`libcudart`, `libcublas`, etc.) directly into the image, so you can compile and run CUDA programs right after pulling them.

But the image alone isn't enough — the driver parts that aren't in the image (`libcuda.so`, the device nodes `/dev/nvidia*`) need to be injected from the host. That requires the NVIDIA Container Toolkit (known as `nvidia-docker2` in the old days), configured as a runtime for Docker:

```bash
# install the toolkit, then configure it as a docker runtime
sudo nvidia-ctk runtime configure --runtime=docker
sudo systemctl restart docker

# run a CUDA container
docker run --rm --gpus all \
  -e NVIDIA_VISIBLE_DEVICES=all \
  nvidia/cuda:12.4.1-base-ubuntu22.04 \
  nvidia-smi
```

### How it works

The core of the old way is the **hook mechanism**. The flow goes roughly like this:

1. `nvidia-container-runtime` is invoked as Docker's runtime; it's essentially a thin wrapper around `runC`;
2. The runtime injects a `prestart` hook into the OCI spec (`config.json`), corresponding to `/usr/share/containers/oci/hooks.d/oci-nvidia-hook.json`;
3. The hook calls `nvidia-container-cli`, which decides which devices to inject based on the `NVIDIA_VISIBLE_DEVICES` environment variable;
4. The CLI mounts/copies the host's device nodes, driver libraries (`libcuda.so`, etc.), and tools like `nvidia-smi` into the container.

So there are a few key points in this mode:

- **Device selection relies on an environment variable**: `NVIDIA_VISIBLE_DEVICES` can be `all`, `0`, `1`, a UUID, and so on;
- **CUDA lives in the image, the driver lives on the host**: the image packages the CUDA user-space libraries, while the driver is injected from the host via the hook; the two are "stitched together" at runtime;
- **Tight coupling to the runtime**: this hook mechanism is NVIDIA-proprietary. Docker needs an `nvidia` runtime configured; Podman recognizes hooks by default, but every runtime has to be adapted separately.

### Pain points of the old way

Over time you run into issues:

- **Huge images**: with the whole CUDA toolchain baked in, images can be several GB, which is a burden for pulling and storage;
- **Version coupling**: the CUDA version in the image and the host driver version need to match; after upgrading the driver, old images may stop working;
- **Too tightly bound to the runtime**: the hook mechanism is NVIDIA's own invention; switching runtimes or to another vendor's devices means starting over;
- **Poor rootless support**: the hook mechanism doesn't play well with rootless containers.

## The new way: CDI (Container Device Interface)

CDI is an open specification under CNCF (`cncf-tags/container-device-interface`) that tries to solve exactly the problems above. Its idea is to standardize "how a device is exposed to a container," so that runtimes, container engines, and orchestrators can all describe devices in the same language.

### Basic concepts

- Devices are identified by **fully-qualified names** in the form `vendor.com/class=unique_name`; for NVIDIA that's things like `nvidia.com/gpu=0`, `nvidia.com/gpu=1:0` (a MIG device), or `nvidia.com/gpu=all`;
- Device information lives in spec files (`.yaml`/`.json`) in the default directories `/etc/cdi` (static configuration) and `/var/run/cdi` (dynamically updated);
- **CDI only makes containers "device-aware"**; resource management (scheduling, quotas) is the orchestrator's job, not CDI's.

A CDI spec file looks roughly like this (excerpt):

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

As you can see, CDI turns "which device nodes to inject, which libraries to mount, which env vars to set" into declarative data. A container engine that sees a device name looks up the spec file and assembles the corresponding content into the OCI spec itself, without going through NVIDIA's proprietary hook.

### How NVIDIA generates CDI specs

The NVIDIA Container Toolkit has supported generating CDI specs since **v1.12.0**. Since **v1.18.0**, a systemd service called `nvidia-cdi-refresh` automatically generates and updates the spec at `/var/run/cdi/nvidia.yaml` on Toolkit install/upgrade, GPU driver install/upgrade, or system reboot:

```bash
# list the currently available CDI devices
nvidia-ctk cdi list
```

One caveat: `nvidia-cdi-refresh` currently does **not** handle driver removal or MIG device reconfiguration automatically. In those scenarios you need to regenerate manually:

```bash
sudo nvidia-ctk cdi generate --output=/var/run/cdi/nvidia.yaml
```

### How to use CDI

**Podman** has native support for specifying CDI devices in the `--device` argument since **v4.1.0**, with no extra configuration:

```bash
podman run --rm \
  --device nvidia.com/gpu=all \
  --security-opt=label=disable \
  ubuntu nvidia-smi -L
```

You can also request a single device by index or by MIG device name:

```bash
podman run --rm \
  --device nvidia.com/gpu=0 \
  --device nvidia.com/gpu=1:0 \
  --security-opt=label=disable \
  ubuntu nvidia-smi -L
```

Note that this runs a plain `ubuntu` image; the CUDA libraries are injected from the host via the CDI spec — which directly addresses the "huge image" pain point of the old way.

**Docker** has supported CDI since **25.0.0** and enables it by default since **28.2.0**. Versions in between (25.0.0 to 28.1.1) require manually enabling it in `/etc/docker/daemon.json`:

```json
{
  "features": {
    "cdi": true
  }
}
```

### What if the runtime doesn't natively support CDI

For runtimes that don't natively support CDI (e.g. Docker with the CDI feature disabled), you can configure the NVIDIA Container Runtime in `cdi` mode. Device selection still uses the `NVIDIA_VISIBLE_DEVICES` environment variable, but the injection goes through CDI:

```bash
sudo nvidia-ctk config --in-place --set nvidia-container-runtime.mode=cdi

docker run --rm -ti --runtime=nvidia \
  -e NVIDIA_VISIBLE_DEVICES=nvidia.com/gpu=all \
  ubuntu nvidia-smi -L
```

One limitation worth remembering: **CDI injection and the NVIDIA Container Runtime hook mechanism are mutually exclusive**. If `/usr/share/containers/oci/hooks.d/oci-nvidia-hook.json` exists, don't use CDI injection at the same time, and don't set `NVIDIA_VISIBLE_DEVICES` either — they'll fight each other.

## Comparing the two approaches

| Dimension | Old way: official CUDA images | New way: CDI |
| --- | --- | --- |
| Device injection | NVIDIA-proprietary prestart hook + `nvidia-container-cli` | Open-spec CDI files parsed by the runtime |
| Device selection | `NVIDIA_VISIBLE_DEVICES` env var | `--device nvidia.com/gpu=xxx` or annotations/CRI fields |
| CUDA source | Baked into the image, huge images | Host driver + CDI injection, plain images work |
| Version coupling | Image CUDA must match driver; upgrade-sensitive | Spec auto-follows the driver via `nvidia-cdi-refresh` |
| Runtime compatibility | Depends on NVIDIA hook, adapted per runtime | Natively unified in Podman/Docker/containerd/CRI-O |
| Rootless | Poor support | Good support |
| Ecosystem openness | NVIDIA-proprietary | CNCF open spec; usable for other vendors' devices too |
| Resource management | No standard | Explicitly left to the orchestrator (e.g. K8s device plugin) |

### Implementation differences

The core difference boils down to one sentence: **the old way is "code injection" (a hook); the new way is "data injection" (a spec file)**.

In the old way, NVIDIA uses a runtime hook to call its own CLI at container start and "shove" devices in; that layer is a black box to users, and every runtime needs separate adaptation. The new way turns device descriptions into standardized data files — the runtime just "reads the spec and assembles the OCI spec," with all the logic in the open ecosystem, so anyone can implement and reuse it.

### Feature and usability differences

- **Image size**: the old way yields multi-GB images; with CDI you can run a plain system image plus the host driver, which is much lighter;
- **Version management**: with the old way, upgrading the driver often means changing the image tag; with CDI, `nvidia-cdi-refresh` automatically follows the driver, so the spec always matches the current driver — far fewer "image/driver mismatch" traps;
- **Device selection**: the old way relies on an obscure environment variable; CDI's fully-qualified names (down to MIG granularity) are intuitive and can go into CRI/annotations, which is much better for expressing intent at the orchestration layer;
- **Multiple runtimes**: one CDI spec works across Podman, Docker, containerd, and CRI-O; with the old way, each runtime may need its own configuration;
- **Rootless**: CDI is friendlier to rootless containers — that was one of its key motivations.

## Summary

The old way is simple and direct. Many people (myself included) ran their first GPU container exactly like that, and it still works today — great for quick validation and personal use. But if your environment is heading toward Kubernetes, needs rootless support, or has multiple runtimes coexisting, CDI is clearly the more modern and less painful direction.

That said, CDI isn't a silver bullet either — `nvidia-cdi-refresh` doesn't cover driver removal or MIG reconfiguration, so the spec sometimes has to be refreshed manually; and it requires runtime-side support, so old Docker versions still need a manual feature flag. But the overall trend is clear: device description is moving from "vendor-proprietary hooks" toward "open, declarative specs."
