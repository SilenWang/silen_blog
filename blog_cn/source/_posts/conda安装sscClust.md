---
title: conda安装sscClust
categories: Others
date: 2018-10-09 08:19:01
tags: ['conda', 'R', 'sscClust', 'Bioconductor', '安装问题']
---

原来又撞到了bug....
<!-- more -->

之前尝试conda安装sscClust失败, 原来是撞到了bug, 在特殊情况下conda安装的R可能会找不到conda安装的bioconductor包, 导致在R里运行时R会重新安装相关依赖.记录完整安装方法如下(假设conda已安装完毕)

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/mro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --set show_channel_urls yes

# 之后进入~/.condarc把里面的'- defualt'删除, 这步非常关键!

conda install -y bioconductor-SingleCellExperiment bioconductor-scran bioconductor-SC3 bioconductor-zinbwave bioconductor-BiocParallel r-base r-devtools r-rcolorbrewer r-rtsne r-class r-factoextra r-cowplot r-data.table r-ggplot2 r-mass r-rjson r-cluster r-ks r-fields r-doparallel r-plyr r-igraph r-densityclust r-e1071
```

之后进入R补充未收录的包

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
source("http://bioconductor.org/biocLite.R")
install.packages(c('RhpcBLASctl', 'ADPclust', 'varSelRF'))
options(unzip = "internal")
devtools::install_github("SilenWang/sscClust", dependencies=F, ref="dev")
```

- 最后devtools::install_github时参数一定要保证不安装依赖, 可能是构建文件有问题, 不管依赖在不在都会强制重新安装依赖?或者是有包的版本不满足?
