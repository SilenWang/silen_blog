---
title: Finally Got the Automated Build Recipe and Action Working with Pixi: Pitfalls and Solutions
categories: Others
date: 2026-02-02 20:51:28
tags: ['pixi', 'recipe', 'github action', 'rattler', 'conda']
---

New tools are always like this: the good parts are amazing, but the unfinished parts can be a headache. Recently, I tried to set up an automated build system using pixi and rattler‑build to regularly package and upload opencode to prefix.dev. The whole process took about 6 hours, during which I encountered quite a few unexpected issues.

<!-- more -->

## Background

I had already successfully created a recipe for opencode and manually uploaded it to prefix.dev. However, opencode updates quite frequently, and I couldn't possibly keep up manually every time. So I wanted to build an automated system similar to bioconda: store the recipe and use GitHub Action to automatically detect new versions, build, and upload conda packages.

It shouldn't be too hard, right? After all, the recipe already exists, writing an Action should be a piece of cake, especially with AI assistance. But the pitfalls I encountered were far more numerous than I had imagined...

## Main Challenges and Solutions

### 1. Cannot use if‑then logic in the recipe's context section

In rattler‑build recipes, the `context` section is used to define variables, but you **cannot** use conditional logic here like you can in the build section. AI also says that the context part can only define variables statically. However, this statement... is both right and wrong. It's true that within YAML syntax, you can only define variables, but... rattler-build actually uses jinja (yes, jinja again) to render a configuration file first, so the supported minijinja functions can implement conditional logic:

```yaml
context:
  architect: ${{ "linux‑x64" if platform == "linux‑64" else "linux‑arm64" }}
```

### 2. Compatibility issues between bump‑recipe and the build section

rattler‑build's `bump‑recipe` feature can automatically update the version and sha256 in the recipe, which is very convenient—I don't need extra steps in the subsequent workflow to handle these operations. However, this feature requires that the recipe's `build` section cannot contain conditional logic. My previous manual recipe used if‑then in the build section, which prevented me from using the bump‑recipe feature.

Therefore, in the updated recipe, I exclusively use minijinja functions and no longer rely on YAML's special syntax.

### 3. AI's incorrect description of environment variable references

When I asked AI how to reference environment variables in a recipe, the answer was to use `${{ VAR_NAME }}`. This is also incorrect; it must be done through the minijinja template engine, with the following syntax:

```yaml
context:
  platform: {{ env.get("TARGET_PLATFORM", default="linux‑64") }}
```

### 4. pixi build can only specify platform; version control requires another approach

I needed to cross‑compile executables for other platforms on a single machine, while also specifying the software version to process. At the same time, I wanted to use `pixi build` instead of `rattler-build build`. This led to another pitfall: `pixi build` is currently an experimental feature, so it cannot pass many parameters to `rattler-build`. For now, the only way to make `pixi` and `rattler-build build` receive consistent parameters is through environment variables:

```bash
TARGET_PLATFORM=linux‑64 VERSION=1.1.35 pixi build
```

Then, in the recipe, retrieve it via `{{ env.get("VERSION") }}`.

### 5. Limitations of variable passing in pixi's built‑in scripts

Because I needed to use `rattler-build bump-recipe`, and to keep the inputs consistent between the action and manual builds (and for simplicity), I set up tasks. Yet another pitfall: pixi can pass arguments to task dependencies, but cannot set environment variables (I needed environment variables to set build parameters). This issue isn't too bad though—it's all shell, so if I can't pass them as environment variables, I can pass them as arguments and then manually set the environment variables:

```toml
[tasks.bump]
cmd = "TARGET_PLATFORM={{ platform }} rattler-build bump-recipe"
args = [    
    { arg = "platform", default = "linux-64" },
]

[tasks.build]
cmd = "TARGET_PLATFORM={{ platform }} pixi build -t {{ platform }}"
args = ["platform"]
depends-on = [{ task = "bump", args = ["{{ platform }}"] }]
```

Note that we need to use `-t` to specify platform again, since `pixi build` seems not using `TARGET_PLATFORM`, or the cross platform build will fail.

### 6. prefix.dev lacks a query API; need to parse pixi search output

Currently, prefix.dev doesn't seem to provide a dedicated API to query version information for specific conda packages. This is easy to solve—just parse the output of the `pixi search` command. Shell as a universal glue is still indispensable...

```bash
version=$(pixi search -q --no‑progress -p linux‑64 -c https://prefix.dev/sylens opencode | grep -oP '\K[0‑9]+\.[0‑9]+\.[0‑9]+' | head -1)
```

### 7. Duplicate downloads with bump‑recipe

`rattler‑build bump‑recipe` downloads the source package once when updating the sha256, and then `pixi build` downloads the same package again. This duplicate download can significantly increase build time when network speed is slow or the package is large. However, there's currently no solution for this, but at least it doesn't affect functionality.


## Complete Implementation

After multiple adjustments, the final pixi project configuration is as follows:

```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda‑forge"]
name = "opencode‑recipe"
platforms = ["linux‑aarch64", "linux‑64"]
version = "0.1.0"
preview = ["pixi‑build"]

[dependencies]
rattler‑build = "*"

[tasks.bump]
cmd = "export TARGET_PLATFORM={{ platform }} && rattler‑build bump‑recipe"
args = [    
    { arg = "platform", default = "linux‑64" },
]

[tasks.build]
cmd = "export TARGET_PLATFORM={{ platform }} && pixi build"
args = ["platform"]
depends‑on = [{ task = "bump", args = ["{{ platform }}"] }]

[package.build.backend]
name = "pixi‑build‑rattler‑build"
version = "0.3.*"
```

Key parts of the recipe file:

```yaml
context:
  name: opencode
  version: {{ env.get("VERSION", default="1.1.35") }}
  platform: {{ env.get("TARGET_PLATFORM", default="linux‑64") }}

package:
  name: {{ name }}
  version: {{ version }}

source:
  url: https://github.com/anomalyco/opencode/releases/download/v{{ version }}/opencode‑{{ "linux‑x64" if platform == "linux‑64" else "linux‑arm64" }}.tar.gz
  sha256: {{ sha256 }}
  file_name: opencode.tar.gz

build:
  number: 0
  string: h{{ hash }}_{{ platform }}
  script: |
    tar xvf opencode.tar.gz
    mkdir -p $PREFIX/bin
    mv opencode $PREFIX/bin/opencode
    chmod 755 $PREFIX/bin/opencode
```

Core steps of the GitHub Action:

```yaml
- name: Run Pixi Task
  run: |
    cd opencode
    
    # Get currently uploaded version
    current_version=$(pixi search -q --no‑progress -p ${{ matrix.platform }} -c https://prefix.dev/sylens opencode | grep -oP '\K[0‑9]+\.[0‑9]+\.[0‑9]+' | head -1)
    
    # Get upstream latest version
    upstream_version="${{ steps.origin_release.outputs.release }}"
    
    if [ "$upstream_version" != "v$current_version" ]; then
      echo "New version detected: $upstream_version (current: $current_version)"
      export VERSION=${upstream_version#v}
      pixi run build ${{ matrix.platform }}
      pixi upload prefix opencode*.conda -c sylens
    else
      echo "Already the latest version"
    fi
```
