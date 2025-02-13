---
title: Writing Dockerfile to Build Docker Image
categories: Script
date: 2019-10-05 16:58:54
tags: ['docker', 'dockerfile']
---

From past experiences... I have gained many skills as a byproduct... Today, this is also the case. My original intention was to use the source code provided by NetEase to compile an improved version of VirtualBox, but in the end, I learned several Docker-related skills...

<!-- more -->

The construction of images is basic to Docker. Although it is possible to manually build images after creating the base image each time, or directly pull them from platforms like Docker Hub, the former is time-consuming and lacks repeatability, while the latter has certain network requirements and limited customizability. Therefore, it is still essential to master the use of Dockerfile, which serves as the blueprint for building images.In fact, this thing is quite simple.

The keywords I might use most frequently in a Dockerfile are:

- `FROM`: Specifies the base image for building. Before starting the build, it will pull the specified image as the base.
- `COPY`: Copies files from the external (host) to the image's internal path. The first parameter is the external path, and the second parameter is the internal path. The usage is similar to `cp`, but I haven't tried using wildcards.
- `ADD`: Similar to `COPY`, but `ADD` can fetch files from the network, which is very useful when deploying containers elsewhere.
- `RUN`: Commands that need to be executed.
- `WORKDIR`: The directory where the built image will start when it runs.
- `CMD`: The command that the built image will run by default when it starts.

Here are some special notes about the `COPY` and `RUN` commands:

- When using `COPY`, you can only fetch files from the build path specified by `docker build` and its subdirectories. This is due to the **context** restriction during Docker runtime.
- Each time the `RUN` command runs, it generates a new layer in the image. If not necessary, try to minimize the number of commands and combine them with `&&` in a single `RUN` command to prevent the image from becoming large and slow.

Here is an example of a Dockerfile I used for deploying an R package `sscClust`:

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

The corresponding R script is:

```r
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages(c("BiocManager", "devtools"))
BiocManager::install(c('SingleCellExperiment', 'scran', 'SC3', 'zinbwave', 'BiocParallel'))
devtools::install_github("Japrin/sscClust")   
```

Execute `docker build -t=test_docker .` to build the image using the above file.

Of course, this example may no longer be able to build successfully because software is constantly updating, and my build process does not specify any specific version of software. As time progresses, various dependency issues may arise... However, as a usage record, it's still quite useful.