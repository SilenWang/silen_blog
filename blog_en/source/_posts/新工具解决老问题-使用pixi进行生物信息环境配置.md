---
title: 新工具解决老问题----使用pixi进行生物信息环境配置
categories: Bioinfomatic
date: 2024-06-22 01:26:08
tags: ['pixi', 'conda']
---

Bioinformatics is an interdisciplinary field, and the toolkit or technology stack used in bioinformatics is also quite "interdisciplinary". The level of fragmentation is, in my opinion, absolutely not less than that of Linux distributions... This also brings us a common challenge: the deployment of bioinformatics analysis environments.

<!-- more -->

## Deploying bioinformatics environments is a hassle

In the development process, one would initially select a programming language, subsequently install the appropriate development environment, and thereafter utilize existing environment or package management tools to handle library dependencies.
However, bioinformatics, like medical statistics I majored in before, was aimed to solve problem first, using whatever is available, and writing codes only if it is needed... Therefore, encountering multiple languages is a highly likely event. Personally, due to the cutting-edge research field of my previous company... I often have to use tools written by university researchers based on the principle of "whatever works, use it", and I have come into contact with a bunch of languages: C, C++, Java, Scala, Perl, Python, R, Go, Julia, Nim, Js... If you also count database languages, configuration file syntax, and the like... The quantity could double again...

Such a messy technology stack makes deploying and setting up a bioinformatics development environment or deploying bioinformatics analysis workflows a very cumbersome task. It's important to know that a significant proportion of tools or software in bioinformatics may no longer be maintained after about a year, which means that the compilers and libraries required by the software will remain outdated, and gradually start to conflict with subsequent versions, causing a more terrifying dependency explosion problem than inter-version upgrades of Linux distributions. Then this dependency is not only going to explode in a single language's toolkit, but there is also a possibility of explosions in toolkits of all languages, and due to the C foundation in Python, R and other scripting languages, there can also be cross explosions... It can make a person's mentality explode in minutes.

However... The problem is far from over, for some adult reasons... Not everyone can get system-level permissions, and when configuring the environment, they can only use the permissions of an ordinary user, so some problems that can be solved by packaging with the distribution's package manager may have to be solved by recompiling the entire software's dependencies... This is enough to make a person explode for a month...

Based on the aforementioned annoying issues, when I first learned about the tool called [conda](https://anaconda.org/), I really felt it was a lifesaver... A large number of bioinformatics software are pre-compiled, packaged, and uploaded, ready to use after download and decompression, and it can automatically calculate if dependencies conflict, ensuring the basic use of tools in different languages, which is simply perfect!

However, after using it for a about a week, you will find that... due to being based on Python, the speed of conda's dependency calculation... is really slow. It might take about 20 minutes to calculate dependencies for around 10 tools, which is especially annoying when there are dependency issues... It's almost like running 20 minutes of calculations, only to find out there's a problem with the variable names, and no breakpoints are set, so you have to start over... Perhaps this downside is so significant that it affects general usage, which led to the development of the C++ implementation of Conda: [Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html), and then came [micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html), and later today's main character, [pixi](https://pixi.sh/latest/).

## Introduction to Pixi

In simple terms, pixi is a multi-language software package manager designed for developers, aimed at making it easy for developers to create and lock down multi-language development environments, making computational results easier to replicate.

The simple usage of pixi (via the command line) can be found in its [official documentation](https://pixi.sh/latest/), and I won't go into detail here because for me, I mainly use it for environment deployment, where a lot of things are installed at once, so I generally use configuration files for this.

## Basic Usage of Pixi

Pixi uses configuration files in the `toml` format (and here we have another one), which is somewhat like `ini`... but with more complex syntax and the ability to support more features.

If you want to manage with `toml`, after installing pixi, you can create a managed project by using `pixi init`. This command will generate a configuration file in the current directory, and the newly generated file will contain the following content.


```toml
[project]
name = "demo"
version = "0.1.0"
description = "Add a short description here"
authors = ["silenwang <silenseek14@gmail.com>"]
channels = ["conda-forge"]
platforms = ["linux-aarch64"]

[tasks]

[dependencies]
```

It already includes the basic information of the project, tasks, and dependencies, which are the most important and main contents. If tasks are not used, they can be left blank. In the dependencies section, you should write down the desired installations (contents will be retrieved from conda). After saving, you can use `pixi install` to install the required contents. During the installation process, a lock file will also be generated. The role of this file is similar to the lock files generated by npm or yarn in nodejs, which is used to record the versions of installed software. By managing this lock file along with the project, you can ensure that the same dependencies are installed again when deploying, guaranteeing the reproducibility of software operation or calculations.

## Handling Multiple Dependencies at Once

Although pixi's software is mainly sourced from conda, it also supports fetching Python packages from pypi. If there are packages not available on conda but can be obtained from pypi, you can add a pypi section in the configuration file. The software within this section will be automatically fetched from pypi and installed after the conda dependencies have been processed, and the obtained packages will also be recorded in the lock file.


```toml
[pypi-dependencies]
spats-shape-seq = "*"
```

If there are additional packages that need to be compiled and installed separately, or further processed with scripts after downloading, you can use the `task` feature to set up tasks.

Tasks set up in the task block can be executed by running `pixi run TASK_NAME`. Although the official documentation does not currently provide a direct way to run specific tasks after dependency processing, you can achieve the goal of automatically proceeding with subsequent compilation and installation tasks after dependency processing by executing the `pixi install -a` command through a task, as shown below:


```toml
[project]
name = "demo"
version = "0.1.0"
description = "Add a short description here"
authors = ["silenwang <silenseek14@gmail.com>"]
channels = ["conda-forge", "bioconda", "milaboratories", "main"]
platforms = ["linux-64"]

[tasks]
install = {cmd = 'pixi install -a'}
deploy = {cmd = 'echo "Deploy done"', depends-on = ["instal"] }

[feature.python3.dependencies]
python = "3.11.*"
pip = "*"
jupyter = "*"
ipykernel = "*"
numpy = "*"
pandas = "*"
statsmodels = "*"
scipy = "*"
tensorflow = "*"
keras = "*"
```

## Deploy Multiple Environments at Once

Pixi also supports the deployment of multiple environments at once, and the general syntax is as follows:


```toml
[feature.python3.dependencies]
python = "3.11.*"
pip = "*"
jupyter = "*"
ipykernel = "*"
numpy = "*"
pandas = "*"
statsmodels = "*"
scipy = "*"
tensorflow = "*"
keras = "*"

[feature.python2.dependencies]
python = "2.7.*"
pip = "*"
jupyter = "*"
ipykernel = "*"

[feature.bioinfo.dependencies]
bwa = "*"
samtools = "*"
bedtools = "*"

[feature.rlang.dependencies]
r-base = "*"
r-irkernel = "*"
r-ggplot2 = "*"
bioconductor-limma = "*"

[environments]
python3 = ["python3"]
python2 = ["python2"]
bioinfo = ["bioinfo"]
rlang = ["rlang"]
```

Specify the environment names in `environments`, and then use blocks like `[feature.rlang.dependencies]` to manage the dependencies within the environment. Once configured, use `pixi install -a` to install all environments (using `pixi install` will only check for dependency issues but will not actually install them).

## Environment Usage

Although pixi uses packages from conda, its environment activation and usage methods follow those of `poetry` and `pipenv`. To activate an environment, run `pixi shell` in the directory; if there are multiple environments, add `-e` to specify the environment name. At the same time, pixi also supports running software within the environment using `pixi run` (yes, `pixi run` is used both for running tasks and for running commands within the environment).
