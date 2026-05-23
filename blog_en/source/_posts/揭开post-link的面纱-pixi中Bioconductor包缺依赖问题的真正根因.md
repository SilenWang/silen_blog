---
title: "Unveiling post-link: The Real Root Cause of Bioconductor Missing Dependencies in pixi"
categories: Bioinformatics
date: 2026-05-22 06:55:00
tags: [pixi, R, Bioconductor, dependency management, conda]
updated: 2026-05-22 06:55:00
---

I previously wrote a blog post about using pixi's tasks feature to fix missing dependency issues with Bioconductor packages (like `GenomeInfoDbData`, `BSgenome.Hsapiens.UCSC.hg38`, etc.) after installation. At the time, I only knew the problem existed but didn't understand the root cause — I was just providing a less-than-ideal workaround.

But recently, with insights from AI responses, I finally figured out the real reason — it all comes down to the **post-link script mechanism** in the Conda ecosystem.

<!-- more -->

## Problem Recap

When managing bioinformatics environments with pixi, after installing R packages that depend on Bioconductor annotation packages (like `Seurat`, `maftools`, etc.), even when these dependencies are explicitly declared in pixi, you still get errors at load time saying the packages don't exist.

A typical error looks like:

```
Error: package 'GenomeInfoDbData' is not installed
```

But in reality, from the generated `pixi.lock`, you can see that `bioconductor-GenomeInfoDbData` and related packages are actually installed.

## The Real Root Cause: post-link Scripts

In the Conda ecosystem, these genome annotation packages are designed to keep package sizes small by **only shipping metadata**. The actual database files (SQLite files, etc.) are downloaded and built locally at install time via **post-link scripts**.

A post-link script is part of the Conda package installation lifecycle — after a package is installed, Conda executes the package's embedded `post-link.sh` (Linux/macOS) or `post-link.bat` (Windows) script to perform post-installation tasks such as:

- Downloading large data files (genome databases, model weights, etc.)
- Compiling native code
- Creating additional directory structures
- Setting environment variables

However, **Pixi, by default, disables (or does not actively execute) Conda package post-link scripts for security and determinism reasons.** This means:

1. Installation completes successfully ✅
2. Package metadata exists ✅
3. But the core data files are never downloaded ❌
4. R naturally can't find the data packages when loading them ❌

This is the real root cause of the "installed but package not found" issue I encountered. And when installing these packages, Pixi actually gives a warning saying that some packages' post-link scripts cannot be executed due to default settings, and you can modify the settings yourself.

Relevant details are documented in [Pixi's documentation](https://pixi.prefix.dev/latest/reference/pixi_configuration/#run-post-link-scripts).

## The Real Solution: Enable post-link Support in Pixi

The official documentation offers two approaches. One is to directly run:

```bash
pixi config set --local run-post-link-scripts
```

to allow these scripts to be executed.

The other approach is essentially the same — manually adding the relevant variable to the configuration file.

This setting only takes effect for packages installed afterward. If the environment has already been configured, you need to run `pixi clean` first, then `pixi install -a` again. You'll notice that after Pixi's progress bar finishes, the log will pause for a while with no further output — that's actually the post-link scripts being executed.

## Comparison with the Previous Approach

The solution mentioned in the earlier blog post {% post_link 修复pixi部署部分bioconda-r包后出现的缺依赖问题 [Approach] %} is still a viable option. However, installing packages directly through R's own dependency system inevitably breaks conda's dependency system, laying hidden traps for future dependency modifications. Therefore, the approach described in this post is the optimal choice.
