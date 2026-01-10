---
title: Making a Small Contribution to Conda-forge
categories: Others
date: 2026-01-09 10:03:56
tags: ['conda-forge', 'tree_sitter_languages', 'PR', 'conda', 'feedstock']
---

Finally, `aider-chat` is available as a conda package, which means it can theoretically be installed globally via `pixi global`. However, during actual installation, you'll find that one of its dependencies, `tree_sitter_languages`, doesn't have a corresponding aarch64 version, causing the installation to fail. This made me wonder: could I rely on AI to solve this?

<!-- more -->

## What is Conda‑forge?

To be honest, before this, I never fully understood the difference between Anaconda and the conda‑forge channel... After reading the documentation, I learned that [Conda‑forge](https://conda-forge.org/) is a community‑driven Conda package distribution channel, while Anaconda is the official resource for conda (and if you use Anaconda casually in a commercial company, you might get warnings or be asked to pay...).
Every package on Conda‑forge corresponds to a **feedstock** repository, which contains all the files needed to build that package, with the most crucial part being the **recipe**.

## What is a recipe?

A recipe is a configuration file that instructs Conda on how to build a software package. It's aptly named a "recipe," and you can find more details in the [official documentation](https://conda-forge.org/docs/maintainer/adding_pkgs).

Since I wasn't building a conda package from scratch this time, I didn't actually need to modify the recipe.

## Preparing the feedstock repository

As mentioned, each software package corresponds to a **feedstock** repository. The [`tree_sitter_languages‑feedstock repository`](https://github.com/conda-forge/tree_sitter_languages-feedstock) is under the conda‑forge organization. According to the official documentation, contributing to a repository is best done by first forking it to your own account, then creating a branch for modifications, and finally pushing changes back to the original repository.

## Modifying the files

Based on Aider's analysis of this feedstock, there's nothing particularly unsupported about aarch64 for this library. It's likely that aarch64 support simply wasn't considered when it was initially set up, and by default, feedstocks don't automatically support aarch64 (probably because there are too few users to justify the effort...). Therefore, the modification is straightforward: just edit `conda-forge.yml` and add:

```yaml
provider:
  linux_aarch64: azure
```

Afterward, following the documentation, I needed to manually run `conda-smithy render -c auto`. Conda would then automatically update the build scripts (truly automated—this is what professional CI/CD looks like!).

## Submitting a PR to the original project for updates

The final step was to commit the code to my forked repository, then submit a Pull Request (PR) from GitHub to the original project. The PR triggers the original repository's CI, which automatically handles building and testing. Once the tests pass and the package maintainer accepts the PR, the package becomes searchable on conda.

With this package, I can now directly run `pixi global install aider-chat` on my Fydetab Duo, and so far, I haven't encountered any issues!
