---
title: docker内进行conda测试
categories: Others
date: 2018-10-06 19:01:12
tags: ['docker', 'conda', 'R', 'sscClust', '容器', '软件安装', '依赖管理', '生物信息学']
---

装软件的时候想了想, 用conda的移植性和可行性还是更好一些, 所以想用conda再尝试一下, 顺便可以看看snakemake怎么样...

<!-- more -->

## 新建docker并安装miniconda

- 初始化一个centos 7的容器, 然后安装`minicoda`, 因为容器是新初始化的所以啥都没...要安装一下必要的东西

```bash
docker run --name "conda_test" -dti IMAGE_ID /bin/bash
docker exec -it DOCKER_ID /bin/bash
yum -y install wget bzip2.x86_64 vim
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc  # conda修改了环境变量, 手动让其生效
```

## 安装R和sscClust

- 本次要模拟的是普通用户使用conda安装sscClust(root获取不可), 因此所有需要的包均不使用yum安装, R包则进入R后进行安装, 不是用conda提供的二进制版本(防止有包conda未收录, 依赖出问题).

- 设置conda的镜像, 然后安装R

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
conda install r
```

- 必要包安装, 由于conda提供的包不包含用于编译的组件, 所以还是要自己编译一部分东西

```bash
wget ftp://xmlsoft.org/libxml2/libxml2-2.7.2.tar.gz ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz https://www.openssl.org/source/openssl-1.1.1.tar.gz https://curl.haxx.se/download/curl-7.61.1.tar.gz

tar xvf openssl-1.1.1.tar.gz
tar xvf gsl-2.5.tar.gz
tar xvf libxml2-2.7.2.tar.gz
```

- 解压后进入相应目录编译安装, 注意编译前加`prefix=/path/to/install`, 然后进行`PATH`配置, 将上述软件安装到的路径加到`PATH`最前面, 这样会优先读取
- curl需要在libgit2之前

- R内安装sscClust

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages("BiocInstaller",repos="https://bioconductor.org/packages/3.7/bioc")
install.packages("devtools")
```

## 测试暂停

- 在测试过程遇到了各种依赖问题, 由于conda提供的包不包含编译需要用的`*-devel`系列包, 如果构建需要的所有包都要从头编译相当花费时间.

## 思路变更

- 对于比较新, 在R源或者bioconductor获取不了的包, 查看其构建的依赖情况, 使用conda安装其依赖项后R内安装该包

```bash

yum install -y unzip wget curl make gcc gcc-gfortran gcc-c++ bzip2.x86_64 vim

conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes

conda install -y -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/ bioconductor-SingleCellExperiment bioconductor-scran bioconductor-SC3 bioconductor-zinbwave bioconductor-BiocParallel

conda install -y r-RColorBrewer r-Rtsne r-class r-factoextra r-cowplot r-data.table r-ggplot2 r-MASS r-rjson r-cluster r-ks r-fields r-doParallel r-plyr r-igraph r-densityClust r-e1071 r-devtools
```

- R内安装的部分

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages(c('varSelRF', 'RhpcBLASctl', 'ADPclust'))
options(unzip = "internal")
source("http://bioconductor.org/biocLite.R")
devtools::install_github('SilenWang/sscClust', dependencies=FALSE, ref="dev")
```

- 使用conda安装的r依然会有依赖问题, 包虽然安装成功, 但是会无法载入(有依赖安装不完全)

## 小结

在无法取得root权限的情况下安装R包会特别麻烦, 还是获取权限装好必要组件后在R内编译安装比较简单, 或者让管理员装了docker之后给docker的使用权限然后自行在docker内操作比较方便.
