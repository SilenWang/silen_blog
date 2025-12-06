---
title: Pixi can now be used in github actions
categories: Coding
date: 2025-12-06 01:43:09
tags: ['github action', 'pixi']
---

In CI/CD workflows, dependency management is often a key factor determining build efficiency and reliability. Recently, I tried the [setup-pixi](https://github.com/prefix-dev/setup-pixi) GitHub Action in a static website deployment pipeline.

<!-- more -->

## Why Choose Pixi?

When I first chose Pixi, it was mainly because it builds on `conda`'s pre‑compiled resources, offers faster dependency resolution than mamba, and has a lock‑file mechanism similar to `nodejs`. These features make it ideal for rapid environment migration and deployment, which are essential in bioinformatics where software changes frequently.

After using it for over a year, I not only manage bioinformatics environments with Pixi but also use it to set up and maintain development environments. Currently, except for system‑resident services like MySQL that it cannot handle, Pixi can install every toolchain I need for my work.

Using a single tool to manage various toolchains, along with its built‑in simple workflow configuration, provides an excellent user experience.

The only minor drawback was a slight inconvenience when I had to use GitHub Actions for automation at my supervisor's request. Since `pixi` still involves virtual environments, installing and running it required some extra environment‑variable tweaks.

However, the official team had already developed the `setup-pixi` [action](https://github.com/marketplace/actions/setup-pixi) long ago—I just didn’t know about it...

## Using `setup-pixi` in GitHub Actions

Using `setup-pixi` is straightforward. Here’s an example:

```yaml
name: DEPLOY

on:
  push:
    branches:
      - main
  workflow_dispatch:

jobs:
  build:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        environment: [website]

    steps:
      - name: Checkout repository
        uses: actions/checkout@v3
        with:
          ref: main

      - name: Setup Pixi
        uses: prefix-dev/setup-pixi@v0.9.3
        with:
          pixi-version: v0.59.0
          environments: ${{ matrix.environment }}

      - name: Run Pixi Task
        run: |
          pixi run build
```

### Step‑by‑Step Explanation

1. **Checkout repository**: Use `actions/checkout` to fetch the latest code. This must come before `Setup Pixi` because the `pixi.toml` and `pixi.lock` files reside in the project.
2. **Setup Pixi**: The `prefix-dev/setup-pixi` Action downloads and installs the specified version of Pixi, then automatically activates the environment defined in `pixi.toml` that matches the one specified in the `matrix` (in this example, the `website` environment).
3. **Run Pixi Task**: With the environment ready, you can freely invoke any task already defined in the project.

The whole process is clean and clear. The Pixi environment ensures dependency consistency during the build phase—define it once, and no further manual setup is needed.

## Summary

By leveraging `setup-pixi` in GitHub Actions, we can easily harness Pixi’s powerful dependency‑management capabilities, simplify configuration, reduce manual setup time, and guarantee high consistency from development to production environments.