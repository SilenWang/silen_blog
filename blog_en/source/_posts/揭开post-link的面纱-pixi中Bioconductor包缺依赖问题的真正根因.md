---
title: Unveiling post-link: The Real Root Cause of Bioconductor Missing Dependencies in pixi
categories: Bioinformatics
date: 2026-05-22 06:55:00
tags: [pixi, R, Bioconductor, dependency management, conda]
updated: 2026-05-22 06:55:00
---

I previously wrote a blog post about using pixi's tasks feature to fix missing dependency issues with Bioconductor packages (like `GenomeInfoDbData`, `BSgenome.Hsapiens.UCSC.hg38`, etc.) after installation. At the time, I only knew the problem existed but didn't understand the root cause — I was just providing a workaround.

Today, with the help of AI-driven deep analysis, I finally figured out the real reason — it all comes down to the **post-link script mechanism** in the Conda ecosystem.

<!-- more -->

## Problem Recap

When managing bioinformatics environments with pixi, after installing R packages that depend on Bioconductor annotation packages (like `Seurat`, `maftools`, etc.), even when these dependencies are explicitly declared in pixi, you still get errors at load time saying the packages don't exist.

A typical error looks like:

```
Error: package 'GenomeInfoDbData' is not installed
```

Yet `conda list` shows that `bioconductor-org.hs.eg.db` and related packages are installed.

## The Real Root Cause: post-link Scripts

In the Conda ecosystem, genome annotation packages like `bioconductor-org.hs.eg.db` are designed to keep package sizes small by **only shipping metadata**. The actual database files (SQLite files, etc.) are downloaded and built locally at install time via **post-link scripts**.

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

This is the real root cause of the "installed but package not found" issue I encountered.

## How to Confirm the Issue

To check if a Bioconductor package relies on post-link for data download, inspect its `extdata` directory:

```bash
# For GenomeInfoDbData as an example
# Find the package installation path
Rscript -e "system.file(package='GenomeInfoDbData')"
# Check the extdata directory — if empty or missing, post-link didn't run
ls $(Rscript -e "cat(system.file(package='GenomeInfoDbData'))")/extdata/
```

If the `extdata` directory doesn't exist or is empty, the post-link script was not executed and the database files were not downloaded.

## Solutions

### Option 1: Use pixi tasks (the previous workaround)

As mentioned in my earlier blog post, use `pixi tasks` to manually install missing data packages via `BiocManager::install()`:

```toml
[tasks]
GenomeInfoDbData = {cmd = 'Rscript -e "BiocManager::install(\"GenomeInfoDbData\")"'}
```

### Option 2: Enable post-link in pixi (recommended)

If you're using pixi >= 0.39.0, you can enable post-link script execution by setting an environment variable:

```bash
export PIXI_ENABLE_POST_LINK=true
```

Or specify it directly when running pixi:

```bash
PIXI_ENABLE_POST_LINK=true pixi install
```

> ⚠️ **Note**: Enabling post-link introduces some security risk, as post-link scripts can execute arbitrary commands. It's relatively safe to use with trusted sources like bioconda, but caution is needed for packages from unknown origins.

### Option 3: Manually execute the post-link script

If you don't want to globally enable post-link, you can run a single package's post-link script manually:

```bash
# Find post-link scripts for the specific package
find $CONDA_PREFIX -name "post-link.sh" -path "*org.hs.eg.db*" -exec bash {} \;
```

### Option 4: Configure pixi.toml to enable post-link

For project-level configuration, you can add to `pixi.toml`:

```toml
[feature]
post-link = true
```

Note that this configuration option's availability varies across pixi versions, so check the documentation for your specific version.

## Summary

| Aspect | Traditional Conda | Pixi (default) | Pixi (post-link enabled) |
|--------|-----------------|----------------|------------------------|
| post-link scripts | Auto-executed | Not executed | Executed |
| Security | Lower | High | Medium |
| Bioconductor data package compatibility | Good | Poor | Good |
| Determinism | Low (network-dependent) | High | Medium |

This investigation shows how important it is to understand the inner workings of your tools. I had been using a workaround for a whole year, essentially working around the problem, when the real solution was right there all along.

If you're using pixi to manage bioinformatics environments and encountering similar Bioconductor missing dependency issues, now you know who the real culprit is — those silent post-link scripts that never got a chance to run.
