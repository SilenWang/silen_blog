---
title: "Pixi Advanced: Using Pixi as a Universal Management Tool for Data Science/Bioinformatics/Web Development"
categories: Coding
date: 2025-12-17 21:56:40
tags: [Pixi, Dependency Management, Data Science, Bioinformatics, Cross‑Platform]
---

Last year, when I joined my current company, I was tasked with preparing a teaching analysis environment for bioinformatics training. This opportunity introduced me to Pixi. Later, as my responsibilities shifted, I not only continued with bioinformatics analysis but also took on full‑stack maintenance work. The projects I handled spanned front‑end and back‑end, involving languages and frameworks beyond the data‑science staples R and Python. Pixi has still been able to serve as a cross‑language/stack development environment management tool (Conda offers an extremely rich resource pool, covering mainstream programming languages and common frameworks). Therefore, I'd like to summarize the practical Pixi features that have proven useful in my daily work.

<!-- more -->

## 1. Handling Shared Dependencies Across Multiple Environments

Currently, about 70‑80% of the bioinformatic projects I work on involve single‑cell analysis. To this day, this field still has two dominant analysis ecosystems: Scanpy in Python and Seurat in R. Both ecosystems continue to produce cutting‑edge new analysis techniques, so using both is unavoidable, and consequently, multiple analysis environments are inevitable. All my projects are highly customized, so installing and testing tools from both ecosystems is a must. Single‑cell data are large, making interactive testing with tools like Jupyter Notebook essential. To avoid repeatedly configuring common packages every time a new environment is created, it's necessary to set Jupyter and the required kernels as shared dependencies.

### Setting a Default Environment

All packages configured under the `[dependencies]` key are placed in the default environment. Every package recorded in this environment will be included in dependency resolution for any other environment.

```toml
[workspace]
authors = ["Sylens Wong <silenseek14@gmail.com>"]
channels = ["conda-forge", "bioconda"]
name = "DEMO1"
platforms = ["linux-64"]
version = "0.1.0"

[environments]
scanpy = ['scanpy']
seurat = ['seurat']

[dependencies]
ipykernel = '*'
r-irkernel = '*'
jupyterlab = '*'

[feature.scanpy.dependencies]
python = '3.*'
scanpy = '*'

[feature.seurat.dependencies]
r-base = '4.*'
r-seurat = '*'
```

### Setting a Common Feature

Another approach is to define a common feature. In `pixi.toml`, we can define multiple environments (`[environments]`), each composed of a set of features. If different environments share common dependencies (e.g., `ipykernel`, `r-irkernel`, `jupyterlab`), we can place them in a feature named `kernel` and then reference that feature in each analysis environment.

```toml
[workspace]
authors = ["Sylens Wong <silenseek14@gmail.com>"]
channels = ["conda-forge", "bioconda"]
name = "DEMO1"
platforms = ["linux-64"]
version = "0.1.0"

[environments]
scanpy = ['kernel', 'scanpy']
seurat = ['kernel', 'seurat']

[feature.kernel.dependencies]
ipykernel = '*'
r-irkernel = '*'
jupyterlab = '*'

[feature.scanpy.dependencies]
python = '3.*'
scanpy = '*'

[feature.seurat.dependencies]
r-base = '4.*'
r-seurat = '*'
```

This way, both the `scanpy` and `seurat` environments automatically include the Jupyter kernel support provided by `kernel`, avoiding duplicate definitions. Moreover, this configuration is more modular: you can configure multiple package components and freely combine them into the environments you need.

## 2. Managing Simple Deployment Tasks and Task Dependencies

Beyond dependency management, Pixi can also define `tasks` and specify dependencies between them. This essentially incorporates a very simple Make‑like workflow system, which is especially useful in projects that require sequential steps such as compilation, packaging, and deployment.

In the `[tasks]` or `[feature.xxx.tasks]` section of `pixi.toml`, we can set up a series of tasks. Each task can declare dependencies using `depends-on`. Pixi ensures that all dependent tasks are executed in order.

Below is an example from a bioinformatics analysis project where several Bioconductor data packages need to be installed sequentially:

```toml
[feature.rplot.tasks]
GenomeInfoDbData = {cmd = 'Rscript -e "BiocManager::install(\"GenomeInfoDbData\")"'}
BSgenome = {cmd = 'Rscript -e "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")"'}
EnsDb = {cmd = 'Rscript -e "BiocManager::install(\"EnsDb.Hsapiens.v86\")"'}
r_dep = {cmd = 'echo "bio dep for R done"', depends-on=['GenomeInfoDbData', 'BSgenome', 'EnsDb']}
```

When you run `pixi run r_dep`, Pixi will automatically execute the three installation tasks first, then output the completion message. Besides this approach of naming tasks as dependencies, Pixi also supports using files as markers to skip already‑completed steps—though the functionality is quite basic, and I haven't used it much myself.

## 3. Specifying Sources for Particular Packages

Sometimes we need to install a specific package from a particular channel (e.g., an internal channel) while other packages still come from the default channels. Pixi allows you to specify a `channel` attribute directly in the dependency declaration. Simply append `{version = "*", channel = "channel-name"}` after the dependency. The only caveat is that `channel-name` must have been defined earlier in the `channels` list.

```toml
[workspace]
name = "DEMO1"
version = "1.1"
channels = ["conda-forge", "bioconda", "dnachun"]
platforms = ["linux-64"]

[feature.label.dependencies]
r-seurat = '5.2.*'
r-SeuratDisk = {version = "*", channel = "dnachun"}
```

Here, `r-SeuratDisk` will be installed from the `dnachun` channel, while other packages follow the default order defined in the `channels` list.

## 4. Mixing PyPI Dependencies

The Conda ecosystem provides a very rich collection of Python resources, but some niche or extremely cutting‑edge packages may still be missing. In such cases, you can additionally specify dependencies from PyPI.

Set the PyPI mirror in `[pypi-options]`, list Conda packages under `[dependencies]`, and list PyPI packages under `[pypi-dependencies]`. If you are using workspace mode, you can also use `[feature.xxx.pypi-dependencies]` within a feature.

```toml
[project]
name = "DEMO1"
version = "1.0"
channels = ["conda-forge", "bioconda"]
platforms = ["linux-64"]

[pypi-options]
index-url = "https://pypi.tuna.tsinghua.edu.cn/simple"

[dependencies]
snakemake = '*'
scanpy = '*'

[pypi-dependencies]
singler = '*'
celldex = '*'
```

In the above configuration, `singler` and `celldex` will be installed via PyPI, while the remaining packages come from Conda.

## 5. Cross‑Platform Dependency Management

In team collaboration or CI/CD scenarios, we may need to run the same environment on Linux‑64, Linux‑aarch64 (my Fydetab Duo uses the Arm‑based RK3588), etc. Pixi can specify different dependency versions per platform, and even use different package sources on different platforms to meet specific requirements (some packages are only available on Conda for Linux‑64).

Use the `[target.<platform>.dependencies]` section to define dependencies for a particular platform. At the same time, you can restrict the platforms an environment supports via the `platforms` field in `[environments]`.

```toml
[workspace]
channels = ["conda-forge"]
platforms = ["linux-64", "linux-aarch64"]

[target.linux-64.dependencies]
git = '*'
git-subtree = '*'

[feature.frontend.dependencies]
nodejs = '18.*'
pnpm = '*'

[feature.frontend.tasks]
front_init = { cmd = "pnpm install", cwd = "app/Frontend_Admin" }
front_build = { cmd = "pnpm run build:prod", depends-on=['front_init']}
```

Here, `git` and `git-subtree` are installed only on `linux-64`; the `frontend` environment is limited to the `linux-64` platform (for example, some front‑end toolchains may not be compatible with ARM).

## 6. Setting Environment Variables for Sub‑environments

Certain tools require environment variables to be set in specific contexts (e.g., Node.js's `NODE_OPTIONS`, proxy settings, etc.). Pixi allows you to define environment variables in a feature via `[feature.xxx.activation.env]`. These variables are automatically set when entering that environment.

```toml
[workspace]
channels = ["conda-forge"]
platforms = ["linux-64", "linux-aarch64"]

[environments]
frontend = ["frontend"]
test = ["test"]

[feature.frontend.activation.env]
NODE_OPTIONS = "--openssl-legacy-provider"

[feature.frontend.dependencies]
nodejs = '18.*'
pnpm = '*'
```

This way, when running tasks inside the `frontend` environment, the `NODE_OPTIONS` environment variable automatically takes effect, avoiding some Node.js version compatibility issues.

## Summary

Through the six typical use cases above, we can see that Pixi is not just a package manager—it's a unified **project environment and task coordination hub**. With a single, concise `pixi.toml` file, it solves dependency management and automation for multi‑language, multi‑platform, multi‑stage tasks, greatly reducing the complexity and maintenance cost of project configuration. Currently, all the projects I write or maintain use Pixi for unified dependency management and task setup. The only limitation is that system‑level services like MySQL cannot be managed through it (but SQLite works fine—the projects aren't that large anyway…).