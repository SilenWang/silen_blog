---
title: 使用singularity部署Rstudio并指定conda环境下的R
tags: ['容器技术', 'Rstudio部署', 'conda环境', 'singularity', 'rstudio-server', 'conda']
categories: Bioinfomatic
date: 2023-07-26 16:13:50
---

今天需要部署个Rstudio给别人用, 记录一下部署的步骤以及要点
<!-- 摘要部分 -->
<!-- more -->

## 前言
Rstudio的运行是必须要靠Root权限来运行的, 不能像Jupyter那样conda或者mamba直装, 因此本次的部署使用了Singularity这种容器(也是我第一次实操使用). Sigularity与Docker或者Podmon在设计上有较大的不同, 主要为HPC环境下的应用而生, 权限管理有较大的不同. 不过这些都不太影响本次使用就是. Sigularity的安装在Rocky Linux下比较简单, 我记得是直接包管理器就能安装了, 因此不作记录

## 1. 容器准备
虽然与Docker是不同的容器技术, Sigularity却提供了完善的Docker镜像兼容支持, 构建好的Docker容器可以很方便的被转换到Sigularity使用的`sif`格式镜像. 本次部署参考的是[rstudio-server-conda](https://github.com/grst/rstudio-server-conda/blob/master/README.md)项目中的说明, 直接使用`singularity pull docker://rocker/rstudio:4.2`拉取docker镜像转换为`sif`镜像使用

## 2. 必要文件准备

使用如下的bash脚本生成需要的若干目录和配置文件

```bash
mkdir -p run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf
echo "auth-minimum-user-id=100" > rserver.conf
echo "session-default-working-dir=/home/silen/Rstudio/Workspace" > rsession.conf
echo "session-default-new-project-dir=/home/silen/Rstudio/Workspace" >> rsession.conf
```

这样一来将生成`run`, `var-lib-rstudio-server`两个文件夹, 以及`database.conf`, `rserver.conf`, `rsession.conf`三个配置文件. 其中`rserver.conf`里的`auth-minimum-user-id=100`配置非常重要, 因为在使用的容器中, Rstudio是以uid为999的rstudio用户来运行的, 如果没有这个用户或者没有切换到该用户的权限, 容器无法正常运行. 因此我参考[网上帖子](https://github.com/grst/rstudio-server-conda/pull/18)中的说法使用`--server-user`参数制定用当前用户运行服务.

然而好巧不巧, 当前用户的uid小于513, 而`auth-minimum-user-id`这一参数默认为1000, 即禁止用户uid小于1000的用户运行服务(似乎是出于安全考虑). 因此就需要在`rserver.conf`里变更参数数值, 来让服务能正常的启动.

## 3. 指定参数启动容器

使用如下命令启动服务即可, 运行的命令主要是将前面创建的目录和配置文件绑定到容器中供服务读取, 另外就是指定了ip与端口方便访问. 由于是担任使用, 就没有弄验证的东西了. 另外还参照[rstudio-server-conda](https://github.com/grst/rstudio-server-conda/blob/master/README.md)项目中的内容制定了使用的conda环境, 这样就不用在容器中重复安装R包了.

```bash
singularity exec \
   --bind run:/run \
   --bind var-lib-rstudio-server:/var/lib/rstudio-server \
   --bind database.conf:/etc/rstudio/database.conf \
   --bind rserver.conf:/etc/rstudio/rserver.conf \
   --bind rsession.conf:/etc/rstudio/rsession.conf \
   --bind /home/silen/R-Plot:/etc/R-Plot \
   --env CONDA_PREFIX=/etc/R-Plot  \
        --env RSTUDIO_WHICH_R=/etc/R-Plot/bin/R \
   rstudio_4.2.sif \
   /usr/lib/rstudio-server/bin/rserver --www-address=0.0.0.0 --www-port=7788 --server-user=$(whoami)
```

启动后界面如下:

![ui](https://raw.githubusercontent.com/SilenWang/Gallary/master/2023/07/upgit_20230727_1690390228.png)

## 4. 注意

Sigularity感觉上会更类似snap和flathub这种软件打包的项目, 宿主机和容器间并不是被那么完全的被隔开. Rstudio启动后, 直接可以访问到宿主机用户的Home目录, 不知道这是这个容器的特殊表现还是Sigularity都是如此? 之后我们需要对内部的分析流程做依赖梳理和容器话... 看来还需要更多的了解Sigularity的原理和使用方式.
