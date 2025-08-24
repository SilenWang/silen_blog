---
title: Build a ChatGLM2 Service Based on the chatglm.cpp Project Using Singularity
tags: ['sinularity', 'container']
categories: Coding
date: 2023-08-16 07:08:46
---

Actually, I've known about singularity for quite some time. As a container specifically designed for HPC (High Performance Computing), it's always been something I wanted to try out. However, just like other NGS technologies outside of Illumina, singularity hasn't gained much traction and it seems even harder to compete with Kubernetes, which is widely adopted by most cloud vendors... Of course, this doesn't really matter to me right now... We're still at least three or five years away from going cloud-based, so using it locally makes perfect sense.

<!-- Abstract section -->
<!-- more -->

So let's start with the most basic step of building an image. The Singularity image build file (Definition File) is quite different from Dockerfile; it's not a script that runs sequentially but rather a configuration file similar to INI files, partitioned into sections for configuration.

The necessary keywords for building an image include:
- `Bootstrap`: Used to define where the base image comes from. There are several types, such as dockerhub or singularity's own library.
- `From`: The name of the base image to build from.
- `%post`: A series of commands used to download files, install software, compile software, etc., during the image build process.

Given that [chatglm.cpp](https://github.com/li-plus/chatglm.cpp) has excellent dependency management, as long as you have a C++ compiler ready, you can install it using pip. Therefore, writing the definition file is quite simple:

```text
Bootstrap: library
From: debian

%post
    apt-get update -y && apt-get upgrade -y
    apt-get install -y cmake gcc g++ pip
    pip install 'chatglm-cpp[api]'
```

Besides these necessary parts, there are other useful sections that might be needed in the future:

- `%files`: Used to copy files into the container during the build process.
- `%environment`: Specifies environment variables for the container.
- `%runscript`: Specifies a script to run when using the `run` subcommand to execute the container.
- `%startscript`: Specifies a script to run when using the `instance` (seems to be used for starting services on boot) subcommand to execute the container.
- `%test`: Specifies commands to run in the container after it's built to check if the build was successful.
- `%labels`: Fills in metadata information for the image.
- `%help`: Provides usage instructions for the image.

With the def file, you can start building by running the command (using fakeroot allows non-root user builds and requires certain settings; see [official documentation](https://docs.sylabs.io/guides/3.11/user-guide/build_a_container.html#fakeroot-builds)). The built image won't be stored locally like Docker images but will generate an sif file.

```bash
singularity build --fakeroot chatglm.cpp.singularity.sif build.def
```

To enable the built image, you can use the following command:

```bash
singularity exec \
   --bind model:/model/ \
   --env MODEL=/model/chatglm2-ggml.q4_0.bin \
   chatglm.cpp.singularity.sif \
   uvicorn chatglm_cpp.openai_api:app --host 0.0.0.0 --port 8877
```
```