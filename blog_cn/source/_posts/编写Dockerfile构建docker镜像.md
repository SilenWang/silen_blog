---
title: 编写Dockerfile构建docker镜像
categories: Script
date: 2019-10-05 16:58:54
tags: ['Docker', 'Dockerfile', '容器化', '镜像构建', 'R语言']
---

从过去的经验来看...我有很多的技能都属于副产物...今天这个也是, 我的本意在于使用网易提供的NemuBox源代码编译出一个网易改进过的Virtualbox, 结果...我学了几样跟docker有关的技能...

<!-- 摘要部分 -->
<!-- more -->

镜像的构建是docker的基本, 虽然也可以每次创建好基础镜像后进入手动构建或者通过类似dockerhub这样的平台直接拉取, 但是前者耗时且重复性差, 后者对网络有一定要求且自定义性差. 掌握编写和修改Dockerfile这一构建图纸的使用方法还是很必要的. 况且...其实这东西也非常简单.

我在Dockerfile中可能会用的比较多的关键字有:

- `FROM`: 指定构建的基础镜像, 在构建开始前会pull指定的镜像作为构建基础.
- `COPY`: 从外部(宿主流经)复制文件到镜像内, 第一个参数是外部路径, 第二个参数是内部路径, 具体参数用法类似`cp`, 似乎不能多个文件同时使用, 没有尝试过通配符是否有效
- `ADD`: 类似`COPY`, 但是`ADD`可以获取网络上的文件, 这在向其他地方部署容器时非常有用
- `RUN`: 需要执行的命令
- `WORKDIR`: 构建好的镜像在启动时进入的目录
- `CMD`: 构建好的镜像在启动时运行的命令

这里`COPY`和`RUN`两个命令需要特别注意:

- `COPY`这个关键字在复制外部路径时, 只能获取`docker build`制定的构建路径下以及其下层的文件. 这个时是受docker运行时**上下文(context)**的限制. 
- `RUN`命令每运行一次都会产生**一层新镜像**, 如果不是必须的, 最好精简需要运行的命令, 并将其以`&&`连接在同一个`RUN`下, 防止构建的镜像又大又慢.

这里我写一个之前工作时需要部署的R包`sscClust`的例子:

```Dockerfile
#
# sscClust Dockerfile
#
# https://github.com/dockerfile/fpm
#

# Pull base image.
FROM centos:7


# Install dep-packages

RUN yum install -y epel-release wget \
    && mv /etc/yum.repos.d/CentOS-Base.repo /etc/yum.repos.d/CentOS-Base.repo.bak \
    && mv /etc/yum.repos.d/epel.repo /etc/yum.repos.d/epel.repo.bak \
    && wget -O /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo \
    && wget -O /etc/yum.repos.d/epel.repo http://mirrors.aliyun.com/repo/epel-7.repo \
    && yum clean all \
    && yum makecache -y \
    && yum install -y R openssl-devel libcurl-devel libxml2-devel gsl-devel \
    && mkdir -p /install_script

# Install R packages

COPY install.r /install_script/install.r

RUN Rscript /install_script/install.r

# Define working directory.
WORKDIR /

# Define default command.
CMD /bin/bash
```

其中涉及的R脚本为:

```r
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c('SingleCellExperiment', 'scran', 'SC3', 'zinbwave', 'BiocParallel'))
devtools::install_github("Japrin/sscClust")   
```

执行`docker build -t=test_docker .`使用上面文件构建镜像

当然现在已经不保证这个例子能顺利构建了, 因为软件总是在不断更新, 而我的构建过程没有任何指定特定版本软件的内容, 随着时间的推进, 很有可能会出现各种依赖爆炸...不过作为使用记录还是阔以哒~
