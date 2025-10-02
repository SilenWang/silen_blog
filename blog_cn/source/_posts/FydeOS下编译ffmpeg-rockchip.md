---
title: FydeOS下编译ffmpeg-rockchip
categories: Others
date: 2025-10-02 16:50:39
tags: ['FydeOS', 'ffmpeg']
---

去年其实就使用Fydetab Duo压制过旅行中的视频，当时找的某个github上给的现成库，结果现在找不到这个库了... 于是这次我又尝试自己来编译了...
<!-- more -->

编译参考了[库作者dnyanmisaka的文档](https://github.com/nyanmisaka/ffmpeg-rockchip/wiki/Compilation)，整个过程还是相对顺利的~

## 使用pixi准备编译环境

这里相比原文档中提到的`git`,`meson`,`cmake`,`pkg-config`,`gcc`,`libdrm-dev`多了很多东西，毕竟fydeos/chromeos并不是传统的linux，缺非常多常见库，因此要补充缺失的依赖程序。

```bash
pixi init
pixi add make gxx meson cmake pkg-config gcc libdrm libdrm-devel-conda-aarch64 pthread-stubs binutils diffutils awk
```

创建pixi环境后，需要`pixi shell`进入虚拟环境，进入后所有后续的编译安装目的地全部指向这个虚拟环境，这样安装后的依赖在最后的ffmpeg编译时才都找的到

## 编译rkmpp
```bash
git clone -b jellyfin-mpp --depth=1 https://github.com/nyanmisaka/mpp.git rkmpp
mkdir -p rkmpp/rkmpp_build
cd rkmpp/rkmpp_build

cmake \
    -DCMAKE_INSTALL_PREFIX=/usr/local/share/project/ffmpeg/.pixi/envs/default/ \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_TEST=OFF \
    ..
make -j $(nproc)
make install
```

## 编译RGA
```bash
git clone -b jellyfin-rga --depth=1 https://github.com/nyanmisaka/rk-mirrors.git rkrga
meson setup rkrga rkrga_build \
    --prefix=/usr/local/share/project/ffmpeg/.pixi/envs/default/ \
    --libdir=lib \
    --buildtype=release \
    --default-library=shared \
    -Dcpp_args=-fpermissive \
    -Dlibdrm=false \
    -Dlibrga_demo=false
meson configure rkrga_build
ninja -C rkrga_build install
```

## 编译ffmpeg并安装
```bash
# ffmpeg部分
git clone --depth=1 https://github.com/nyanmisaka/ffmpeg-rockchip.git ffmpeg
cd ffmpeg
./configure --prefix=/usr/local/share/project/ffmpeg/.pixi/envs/default \
    --enable-gpl --enable-version3 --enable-libdrm --enable-rkmpp --enable-rkrga
make -j $(nproc)
make install
```

## 后记
执行上面所有的步骤后，虚拟环境内就有rockchip芯片专用的ffmpeg了，视频处理速度比cpu版本的快了10x有多~

晚点考虑把上面的内容做一个带任务的Pixi环境，然后也学习下怎么用pixi打包和上传包。