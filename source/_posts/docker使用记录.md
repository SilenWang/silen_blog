---
title: docker使用记录
categories: Others
date: 2018-09-02 20:15:02
tags: ['docker']
---

为了充分利用二手不出去的烂设备...我决定自己搭服务器积累经验...考虑到搭建的东西的移植性, 使用docker放在容器内执行, 这样以后方面迁移与使用.

<!-- more -->

# 安装

环境系统使用manjaro, 毕竟这个我最熟悉. 并且archwiki能为许多奇怪的问题提供帮助, 安装docker的话很简单:

```bash
sudo pacman -S docker
# 为了让普通用户使用docker, 需要添加相应用户组并把要使用的用户加入组中
sudo groupadd docker
sudo gpasswd -a ${USER} docker
# 上面用usermod命令也可以:
sudo usermod -Ga docker ${USER}
# 如果已经启用过则`restart`, 不需要`enable`
systemctl start docker
systemctl enable docker
```

为加快镜像(images)拉取速度, 设施dockerhub的国内镜像源(registry mirror), 配置文件是`/etc/docker/daemon.json`, 配置完成后执行`sudo systemctl restart docker`:

```txt
{
  "registry-mirrors": ["https://docker.mirrors.ustc.edu.cn/"]
}
```

# 镜像创建及镜像内部操作

准备在镜像内安装R并安装`sscClust`包, 公司的服务器是CentOS6, 所以使用拉取相应镜像:

```bash
docker pull centos:6
```

但是在进行了上述拉取之后, 一直不能正常的运行镜像, 根据[这里](https://forums.docker.com/t/docker-run-it-has-started-failing-with-status-139/18309)的描述, 似乎这是官方镜像有问题, 所以改为拉取了`centos:7`, 确实在更换版本后正常启动.

```bash
docker run -dti IMAGE_ID /bin/bash
# 登陆入dcoker进行操作, 注意, 登陆后是root, 操作小心...
docker exec -it CONTAINER_ID /bin/bash
```

登陆后安装vim, 对环境变量作基本设置, 然后开始安装R

```bash
yum -y install epel-release
# 更改镜像源
yum install wget
mv /etc/yum.repos.d/CentOS-Base.repo /etc/yum.repos.d/CentOS-Base.repo.bak
mv /etc/yum.repos.d/epel.repo /etc/yum.repos.d/epel.repo.bak
wget -O /etc/yum.repos.d/CentOS-Base.repo http://mirrors.aliyun.com/repo/Centos-7.repo
wget -O /etc/yum.repos.d/epel.repo http://mirrors.aliyun.com/repo/epel-7.repo
# 不同的镜像源包版本会不一样, 删除已有缓存防止报错
yum clean all
yum -y makecache
yum -y install R
```

补充R包安装需要的一些系统软件包

```bash
mkdir /usr/share/doc/R-3.5.0/html
yum install -y openssl-devel libcurl-devel libxml2-devel gsl-devel
```

之后进入R开始安装`sscClust`(安装R本来就会自动安装一大堆依赖, 安装这个包的时候会有更多依赖...)

```r
options(repos="https://mirrors.shu.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.ustc.edu.cn/bioc/")
source("http://bioconductor.org/biocLite.R")
install.packages("devtools")
biocLite(c('SingleCellExperiment', 'scran', 'SC3', 'zinbwave', 'BiocParallel'))
devtools::install_github("Japrin/sscClust")
```

如果网络没有特别问题, `sscClust`应该是能够顺利安装的, 最后将docker容器的更改保存为一个镜像(image), 保存完成后会显示一串sha256值, 使用`docker images`可看到已经保存下来的镜像

```bash
docker commit CONTAINER_ID r/sscclust
```