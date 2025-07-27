---
title: 使用singularity构建一个基于chatglm.cpp项目的ChatGLM2服务
tags: ['容器技术', 'singularity', 'ChatGLM2', 'chatglm.cpp', '容器构建', 'HPC']
categories: Bioinfomatic
date: 2023-08-16 07:08:46
---


其实很早之前就知道singularity这个东西了，作为有别于docker，专门为HPC开发的容器，一直都想试试。奈何就像Illumina以外的其他NGS技术一样，singularity虽然没有挂，但是至今都没有什么声量，而且在k8s这种容器集群管理方案被绝大多数云厂商采纳后，singularity要竞争似乎更难了... 当然这跟我现在没什么关系... 咱目前离上云感觉至少还有三五年的距离，因此在本地集群启用这个感觉完全合乎情理。

<!-- 摘要部分 -->
<!-- more -->

那么还是从最基础的构建镜像开始，singularity的镜像构建文件（Definition File）和Dockerfile还是很不一样的，不是另一种从前往后执行的脚本，感觉更像是ini这种分区进行配置的配置文件。

构建一个镜像所必要的关键字包括：
- `Bootstrap`: 用于定义base镜像从哪里来，有多种类型，可以是dockerhub可以是singularity维护的镜像lirary
- `From`: 构建base镜像的名称
- `%post`: 构建镜像时下载文件，安装软件，编译软件等的一系列命令

鉴于[chatglm.cpp](https://github.com/li-plus/chatglm.cpp)本身的依赖做的非常好，其实只要准备好c++编译器，用pip命令就能装了，所以写起来特别简单：

```text
Bootstrap: library
From: debian

%post
    apt-get update -y && apt-get upgrade -y
    apt-get install -y cmake gcc g++ pip
    pip install 'chatglm-cpp[api]'
```

除了这些必要的内容，以后可能还用得上的东西包括：

- `%files`：构建时向容器内拷贝文件用的部分
- `%environment`：用于规定换纪念馆变量
- `%runscript`：用于指定使用`run`子命令运行容器时会使用的脚本
- `%startscript`：用于指定使用`instance`（似乎是开机启动服务）子命令运行容器时会使用的脚本
- `%test`：指定容器构建完成后在容器内运行的命令（用于检查构建是否完成）
- `%labels`：填写镜像的meta信息
- `%help`：填写镜像的使用说明


有了def文件，接下来运行命令就可以开始构建（使用fakeroot是使用非root用户构建，还需要作一定设置参见[官方文档](https://docs.sylabs.io/guides/3.11/user-guide/build_a_container.html#fakeroot-builds)），构建完成的镜像不像docker是自己存储起来，而是会生成一个sif文件。

```bash
singularity build --fakeroot chatglm.cpp.singularity.sif build.def
```

构建好的镜像用如下命令就可以启用啦：

```bash
singularity exec \
   --bind model:/model/ \
   --env MODEL=/model/chatglm2-ggml.q4_0.bin \
   chatglm.cpp.singularity.sif \
   uvicorn chatglm_cpp.openai_api:app --host 0.0.0.0 --port 8877
```
