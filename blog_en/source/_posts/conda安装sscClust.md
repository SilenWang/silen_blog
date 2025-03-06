---
title: conda Install sscClust
categories: Script
date: 2018-10-09 08:19:01
tags: ['conda']
---

It seems I hit a bug again...
<!-- more -->

I previously tried to install `sscClust` using conda, but it failed because of a bug. In certain cases, the R installed via conda might not be able to find the bioconductor packages installed via conda, leading to R reinstalling the dependencies when running. I recorded the complete installation method as follows (assuming conda is already installed):

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/mro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --set show_channel_urls yes

# After that, delete the '- defualt' in ~/.condarc. This step is very crucial!

conda install -y bioconductor-SingleCellExperiment bioconductor-scran bioconductor-SC3 bioconductor-zinbwave bioconductor-BiocParallel r-base r-devtools r-rcolorbrewer r-rtsne r-class r-factoextra r-cowplot r-data.table r-ggplot2 r-mass r-rjson r-cluster r-ks r-fields r-doparallel r-plyr r-igraph r-densityclust r-e1071
```

After that, supplement the packages not included in R:

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
source("http://bioconductor.org/biocLite.R")
install.packages(c('RhpcBLASctl', 'ADPclust', 'varSelRF'))
options(unzip = "internal")
devtools::install_github("SilenWang/sscClust", dependencies=F, ref="dev")
```

- Finally, when using `devtools::install_github`, make sure to set the parameter to not install dependencies. It might be due to issues with the build files, forcing a reinstallation of dependencies? Or is it because the package versions do not meet the requirements?
