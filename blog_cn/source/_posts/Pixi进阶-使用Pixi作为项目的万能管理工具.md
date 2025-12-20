---
title: Pixi进阶：使用Pixi作为数据科学/生物信息学/web开发的万能管理工具
categories: Coding
date: 2025-12-17 21:56:40
tags: [Pixi, 依赖管理, 数据科学, 生物信息学, 跨平台]
---

去年我入职新公司的时候，接到了一项为生物信息培训准备教学用分析环境的任务，以此为契机，我接触到了Pixi。后来由于工作内容变动，我不仅要进行生物信息分析，还要兼任全栈维护的工作，接手的项目横跨前段后端，涉及的语言和框架不再限于数据科学常用的R、Python，Pixi也依然能充当一个跨语言/技术栈的开发环境管理工具（Conda资源极丰富，主流的编程语言和常用框架都有资源），因此现在再次来小结下，在我的工作中用得上的Pixi实用功能。

<!-- more -->

## 1. 处理多环境的公用依赖

现在经手的项目七八成都是单细胞，这个领域到目前为止，都还是有两个主流的分析生态，即Python下的Scanpy，以及R下的Seurat，并且两种生态下都在持续产出前沿的新分析技术，因此不可避免的两个生态都要用，分析环境也就必然是多个。而我做的项目全部都是个性化的，因此安装和测试两套生态下的工具不可避免，单细胞数据又大，使用 Jupyter Notebook 这种工具进行交互性测试是必然的，为了避免每次创建新环境都反复配置共用软件包，将 Jupyter 以及必要的内核设为通用依赖是很必要的。

### 设置defualt环境

在设置环境的关键字`[dependencies]`下配置的所有软件包都在defualt环境中，这个环境中记录的所有包都会在任一环境中被纳入依赖计算。

```toml
[workspace]
authors = ["Sylens Wong <qium@aimingmed.com>"]
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

### 设置通用feature

另一种方式，则是设置通用feature。在`pixi.toml`中，我们可以定义多个环境（`[environments]`），每个环境由一组feature组成。如果不同的环境有公共依赖（如`ipykernel`、`r-irkernel`、`jupyterlab`）放入一个名为`kernel`的feature中，然后在各个分析环境中引用该feature。

```toml
[workspace]
authors = ["Sylens Wong <qium@aimingmed.com>"]
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

这样，无论是`scanpy`环境还是`seurat`环境，都会自动包含`kernel`提供的Jupyter内核支持，避免了重复定义。另外，这种配置方式更模块化，可以配置多个包组件，然后自由组合成需要的环境。

## 2. 管理简单的部署任务以及任务之间的依赖

除了依赖管理，Pixi还可以定义任务（`tasks`），并指定任务之间的依赖关系。等于是纳入了一个非常简单的Make工作流体系。这在需要顺序执行编译、打包、部署等步骤的项目中尤其有用。

在`pixi.toml`的`[tasks]`或`[feature.xxx.tasks]`节中，我们可以设置一系列的任务，每一个任务可以利用`depends-on`声明任务依赖。Pixi会确保所有依赖任务按顺序执行。

下面是一个生物信息分析项目中，需要依次安装几个Bioconductor数据包的任务定义：

```toml
[feature.rplot.tasks]
GenomeInfoDbData = {cmd = 'Rscript -e "BiocManager::install(\"GenomeInfoDbData\")"'}
BSgenome = {cmd = 'Rscript -e "BiocManager::install(\"BSgenome.Hsapiens.UCSC.hg38\")"'}
EnsDb = {cmd = 'Rscript -e "BiocManager::install(\"EnsDb.Hsapiens.v86\")"'}
r_dep = {cmd = 'echo "bio dep for R done"', depends-on=['GenomeInfoDbData', 'BSgenome', 'EnsDb']}
```

执行`pixi run r_dep`时，Pixi会自动先执行三个安装任务，最后输出完成信息。除了这种指定任务名为依赖的组织方式，Pixi也支持以文件作为标志，从而跳过已完成步骤的组织方式，不过... 功能过于简单，我自己没有太多使用。

## 3. 为特定软件包指定来源

有时我们需要从特定的channel安装某个软件包（例如某个内部channel），而其他包仍从默认channel获取。Pixi允许在依赖声明中单独指定`channel`属性。在依赖项后面添加`{version = "*", channel = "channel-name"}`即可。唯一要注意的是，`channel-name`必须在前面的`channel`中被定义过。

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

这里`r-SeuratDisk`将从`dnachun`这个channel安装，而其他包仍使用`channels`列表中定义的默认顺序。

## 4. 混合使用PyPI依赖

Conda生态提供了非常丰富的Python资源，但是一些冷门或特别前沿的包还是会找不到的，这时候可以额外指定Pypi来源的依赖。

在`[pypi-options]`中设置PyPI镜像，然后在`[dependencies]`中列出conda包，在`[pypi-dependencies]`中列出PyPI包。如果使用workspace模式，也可以在feature中使用`[feature.xxx.pypi-dependencies]`。

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

上述配置中，`singler`和`celldex`将通过PyPI安装，其余包通过conda安装。

## 5. 跨平台依赖管理

在团队协作或CI/CD中，我们可能需要在Linux‑64、Linux‑aarch64（我的Fydetab Duo是Arm架构的RK3588）运行同一套环境。Pixi可以针对不同平台指定不同的依赖版本，甚至在不同的平台使用不同的包来源，以满足特定的依赖（有些包Conda上只有Linux64）。

使用`[target.<platform>.dependencies]`节可以为特定平台定义依赖。同时，可以在`[environments]`中利用`platforms`字段限制环境支持的平台。

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

这里，`git`和`git-subtree`只在`linux-64`上安装；而`frontend`环境仅支持`linux-64`平台（例如某些前端工具链在ARM上可能不兼容）。

## 6. 子环境的环境变量设定

某些工具需要在特定环境下设置环境变量（例如Node.js的`NODE_OPTIONS`、代理设置等）。Pixi允许在feature中通过`[feature.xxx.activation.env]`定义环境变量，这些变量会在进入该环境时自动设置。

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

这样，当运行`frontend`环境中的任务时，`NODE_OPTIONS`环境变量会自动生效，避免一些Node.js版本兼容性问题。

## 总结

通过以上六个典型案例，我们可以看到Pixi不仅仅是一个包管理器，更是一个统一的**项目环境与任务协调中心**。它用一份简洁的`pixi.toml`文件，解决了多语言、多平台、多阶段任务的依赖管理与自动化问题，大大降低了项目配置的复杂度和维护成本。目前我所有自己编写和接收维护的项目都用它来统一进行依赖管理和任务设置。唯一的限制是MySQL这种一定要系统级服务的工具不可能通过它来管理（不过用SQLite就是了，项目其实都不大...）