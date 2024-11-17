---
title: 新工具解决老问题----使用pixi进行生物信息环境配置
categories: Bioinfomatic
date: 2024-06-22 01:26:08
tags: ['pixi', 'conda']
---

生物信息是一门交叉学科，而生物信息学中使用的工具集或技术栈，也是相当的“交叉”，其片化的程度，我是觉得绝对不亚于Linux发行版... 这也带来了大家都会面临的难题：生物信息分析环境部署。

<!-- 摘要部分 -->
<!-- more -->

## 生信环境部署是个麻烦事

一般来说，用什么语言做开发，那就部署那个语言的开发环境，然后使用相应的环境/包管理软件对库的依赖管理就好了。但是生物信息跟我之前学的医学统计一样，应用成分重，有啥用啥，没有再写... 因此，接触多种语言这种事是极大概率事件。就我个人来说，由于前公司的研究领域教前沿... 很多时候都要使用大学里的研究者们基于`会啥用啥`这一抓到老鼠就是好猫一般的原则下写出来的工具，我已经接触过一堆语言了：C、C++、Java、Scala、Perl、Python、R、Go、Julia、Nim、Js... 如果再要算上一些数据库语言、配置文件语法之类的... 数量还能再翻个一倍...

这么杂乱的技术栈，就让部署和设置生物信息开发环境，或者部署生物信息分析流程成了一件非常繁琐的事情。要知道生物信息学中有相当比例的工具或软件，在一年左右就不一定会再有人去维护了，这意味着软件需要的编译器和库会一直停留在过时的版本，然后逐渐开始和后续版本打架，造成比Linux发行版跨版本升级更可怕的依赖爆炸问题。然后这个依赖还不是只有一个语言下的工具集会爆炸，是所有语言下的工具集都有爆炸的可能，且由于Python、R等一些脚本语言底层有C，所以还能交叉爆炸... 分分钟搞得人心态爆炸。

然而... 为题远不止步于此，处于一些大人的原因... 不是所有人都能拿到系统级别的权限，配置环境时仅能使用普通用户的权限，因此一些用发行版的包管理器装个包就能解决的问题，最后可能只能选择将整个软件的依赖从头编译一遍... 这又够人爆炸一个月的了...

基于前述种种令人恼火的问题，当我第一次知道有[conda](https://anaconda.org/)这么个工具的时候，真的是感觉救我狗命啊... 大量的生物信息学软件都已预编译后打包上传，下载解压缩即用，还能自动计算依赖是否打架，保证不同语言下工具的基本使用，简直是太完美了！

不过用过一周之后就会了解到... 由于基于python，conda的依赖计算速度... 是真的慢，10个左右的工具可能就会算个20分钟了，这在依赖会有问题的时候尤其令人恼火... 约等于跑了20分钟的计算，发现变量名称写得有问题，又没有设好断点，只能从头再来一遍... 可能也就是这个缺点大到影响一般使用，就有人开发了Conda的C++实现：[Mamba](https://mamba.readthedocs.io/en/latest/user_guide/mamba.html)，然后又有了[micromamba](https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html)，再后来有了今天的主角，[pixi](https://pixi.sh/latest/)

## pixi简介

简单来说，pixi是一个为开发者设计的多语言软件包管理程序，目的在于使开发者能够轻松的创建和锁定多语言开发环境，使计算结果更容易被重复。

pixi的简单使用方法(通过命令行)可以参见它的[官方文档](https://pixi.sh/latest/)，这里不多说，因为对我来说主要会用它进行环境部署，会一次性安装一大堆的东西，因此一般都是用配置文件来进行配置的。

## pixi基本使用

pixi使用`toml`格式的配置文件(得，又多了一种)，其格式有点类似ini...但是语法更复杂，能支持更多的东西。

想使用`toml`来管理的话，安装好pixi后，使用`pixi init`就可以创建一个管理项目，并在当前目录生成一份配置文件了，刚生成的文件包含如下内容。

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

其中已经包含了项目基本信息、任务、依赖这三块最重要，也最主要的内容。任务不使用直接留空就可以了，在依赖块内写上需要的安装的内容（内容将从conda获取），保存后就可以使用`pixi install` 来安装需要的内容了，安装过程中还会生成一个lock文件，这个文件的作用跟nodejs中npm或yarn生成的lock文件作用类似，用于记录被安装的软件版本。将这个lock文件随项目一起管理，就能保证再次部署时能安装一模一样的依赖，保证软件可运行或计算可重复了。

## 多种依赖一次性处理

pixi的软件虽然主要从conda获取，但是它同时支持从pypi获取python包。如果有conda上没有，但是可以从pypi获取的包，可额外在配置文件中增加pypi区块。区块内的软件将会在conda依赖处理完成后自动从pypi获取软件包并进行安装，获取的包同样会记录到lock文件中。

```toml
[pypi-dependencies]
spats-shape-seq = "*"
```

如果还有其他包要单独进行编译安装，或者下载后用脚本进一步处理，则可以使用`task`特性来设置任务。

在task块中设置好的任务可以通过`pixi run TASK_NAME`的方式来进行运行。当前官方虽然没有给在依赖处理完成后，直接运行特定任务的方式，但是可以通过下面的写法来通过任务执行`pixi install -a`命令以先处理依赖，从而达到在依赖处理完成后自动进行后续编译安装任务的目的：

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

## 多种环境一次性部署

pixi还支持了多环境一次性部署，大概的写法如下：

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

即在`environments`中指定环境的名称，然后使用名字类似`[feature.rlang.dependencies]`的块来处理环境内的依赖情况。配置好后，使用`pixi install -a`来安装所有的环境（使用`pixi intall`只会检查依赖是否有问题，但是不会实际安装）

## 环境使用

pixi虽然使用来自conda的包，但是其环境激活和使用方式却采用了`poetry`和`pipenv`的方式，激活环境需要在目录下运行`pixi shell`，如果有多个环境，则加`-e`指定环境名。同时pixi也支持`pixi run`的方式来运行环境内的软件（对，`pixi run`既用来运行任务，也用来运行环境内命令）。

## 给个栗子

施工中...