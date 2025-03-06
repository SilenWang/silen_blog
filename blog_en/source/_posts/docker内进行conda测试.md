---
title: Testing Conda Inside Docker
categories: Script
date: 2018-10-06 19:01:12
tags: ['docker', 'conda']
---

When installing software, I thought that the portability and feasibility of using conda were better, so I decided to try it out and see how Snakemake works...

<!-- more -->

## Creating a Docker Container and Installing Miniconda

- Initialize a CentOS 7 container and install `minicoda` because the container is newly initialized, so there's nothing installed... need to install some necessary things.

```bash
docker run --name "conda_test" -dti IMAGE_ID /bin/bash
docker exec -it DOCKER_ID /bin/bash
yum -y install wget bzip2.x86_64 vim
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
sh Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc  # conda modified the environment variables, manually make it effective
```

## Installing R and sscClust

- This time, we want to simulate a regular user installing `sscClust` (root access is not available), so all necessary packages will not be installed using yum. R packages will be installed within R itself, not using the binary versions provided by conda (to prevent issues with packages not included in conda or dependency problems).

- Set up the conda mirror and then install R.

```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
conda config --set show_channel_urls yes
conda install r
```

- Install necessary packages. Since conda does not provide the compilation components, we still need to compile some parts ourselves.

```bash
wget ftp://xmlsoft.org/libxml2/libxml2-2.7.2.tar.gz ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz https://www.openssl.org/source/openssl-1.1.1.tar.gz https://curl.haxx.se/download/curl-7.61.1.tar.gz

tar xvf openssl-1.1.1.tar.gz
tar xvf gsl-2.5.tar.gz
tar xvf libxml2-2.7.2.tar.gz
```

- Extract and enter the corresponding directories to compile and install, note that before compiling, add `prefix=/path/to/install`, then configure the `PATH` to add the installation path of the above software at the beginning of `PATH`, so it will prioritize reading from there.
- curl needs to be installed before libgit2.

- Install sscClust within R.

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages("BiocInstaller",repos="https://bioconductor.org/packages/3.7/bioc")
install.packages("devtools")
```

## Testing Pause

- During the testing process, various dependency issues were encountered. Since conda does not provide the compilation components needed for `*-devel` series packages, if all build dependencies need to be compiled from scratch, it will take a lot of time.

## Change in Approach

- For newer packages that cannot be obtained from R source or Bioconductor, check their dependency situation and install the dependency items using conda before installing the package within R.

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

- Install some packages within R.

```R
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
install.packages(c('varSelRF', 'RhpcBLASctl', 'ADPclust'))
options(unzip = "internal")
source("http://bioconductor.org/biocLite.R")
devtools::install_github('SilenWang/sscClust', dependencies=FALSE, ref="dev")
```

- Even though the R package installed using conda still has dependency issues, it will be unable to load (dependencies are not fully installed).

## Conclusion

Installing R packages without root access is particularly troublesome. It's simpler to install necessary components first and then compile within R itself, or for administrators to set up Docker and give usage permissions to the user, allowing them to operate within Docker themselves.
```