---
title: Docker Usage Record
categories: Others
date: 2018-09-02 20:15:02
tags: ['docker']
---

To make good use of the old, broken device... I decided to set up a server to accumulate experience. Considering the portability of the things I will be setting up, using Docker in containers is convenient for migration and usage.

<!-- more -->

# Installation

The environment system uses Manjaro, since this is what I am most familiar with. And ArchWiki can provide help for many strange problems. Installing Docker is very simple:

```bash
sudo pacman -S docker
# To allow ordinary users to use Docker, you need to add the corresponding user group and add the user to the group.
sudo groupadd docker
sudo gpasswd -a ${USER} docker
# Alternatively, you can use the `usermod` command:
sudo usermod -Ga docker ${USER}
# If it has already been enabled, restart it; no need to enable it again.
systemctl start docker
systemctl enable docker
```

To speed up image (images) downloads, set up a domestic mirror source (registry mirror) for Docker Hub. The configuration file is `/etc/docker/daemon.json`. After configuring, execute `sudo systemctl restart docker`:

```txt
{
  "registry-mirrors": ["https://docker.mirrors.ustc.edu.cn/"]
}
```

# Image Creation and Internal Operations

Prepare to install R and the `sscClust` package in a container. The company's server is CentOS6, so pull the corresponding image:

```bash
docker pull centos:6
```

However, after pulling the above image, it could not start normally. According to [here](https://forums.docker.com/t/docker-run-it-has-started-failing-with-status-139/18309), this might be an issue with the official image. So I switched to pulling `centos:7`, which worked fine after changing versions.

```bash
docker run -dti IMAGE_ID /bin/bash
# Log in to Docker for operations. Note that when logged in, you are root; operate carefully...
docker exec -it CONTAINER_ID /bin/bash
```

After logging in, install vim and set up basic environment variables. Then start installing R.

```bash
yum -y install epel-release
# Change the image source
yum install wget
mv /etc/yum.repos.d/CentOS-Base.repo /etc/yum.repos.d/CentOS-Base.repo.bak
mv /etc/yum.repos.d/epel.repo /etc/yum.repos.d/epel.repo.bak
wget -O /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo
wget -O /etc/yum.repos.d/epel.repo http://mirrors.aliyun.com/repo/epel-7.repo
# Different image sources have different package versions, delete existing cache to prevent errors.
yum clean all
yum -y makecache
yum -y install R
```

Install some system software packages required for R package installation.

```bash
yum install -y openssl-devel libcurl-devel libxml2-devel gsl-devel
```

After entering R, start installing `sscClust` (installing R will automatically install a lot of dependencies; installing this package will have even more dependencies...)

```r
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
source("http://bioconductor.org/biocLite.R")
install.packages("devtools")
biocLite(c('SingleCellExperiment', 'scran', 'SC3', 'zinbwave', 'BiocParallel'))
devtools::install_github("Japrin/sscClust")
```

If there are no special network issues, `sscClust` should install smoothly. Finally, save the changes made to the Docker container as an image. After saving, it will display a series of sha256 values. Use `docker images` to see the saved images.

```bash
docker commit CONTAINER_ID r/sscclust
```
```