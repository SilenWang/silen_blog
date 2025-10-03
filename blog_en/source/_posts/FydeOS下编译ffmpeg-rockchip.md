---
title: Compiling ffmpeg-rockchip on FydeOS
categories: Others
date: 2025-10-02 16:50:39
tags: ['FydeOS', 'ffmpeg']
---

Last year, I actually used Fydetab Duo to compress videos from my travels. At that time, I found a pre-built library on GitHub, but now I can't find that repository anymore... So this time, I tried to compile it myself...

<!-- more -->

The compilation referred to [the documentation by library author dnyanmisaka](https://github.com/nyanmisaka/ffmpeg-rockchip/wiki/Compilation), and the entire process was relatively smooth~

## Preparing the compilation environment using pixi

Compared to the `git`, `meson`, `cmake`, `pkg-config`, `gcc`, `libdrm-dev` mentioned in the original documentation, many more things are needed here. After all, FydeOS/ChromeOS is not a traditional Linux distribution and lacks many common libraries, so we need to supplement the missing dependencies.

```bash
pixi init
pixi add make gxx meson cmake pkg-config gcc libdrm libdrm-devel-conda-aarch64 pthread-stubs binutils diffutils awk
```

After creating the pixi environment, you need to enter the virtual environment with `pixi shell`. Once inside, all subsequent compilation and installation destinations will point to this virtual environment, ensuring that all dependencies installed can be found during the final ffmpeg compilation.

## Compiling rkmpp
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

## Compiling RGA
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

## Compiling and installing ffmpeg
```bash
# ffmpeg part
git clone --depth=1 https://github.com/nyanmisaka/ffmpeg-rockchip.git ffmpeg
cd ffmpeg
./configure --prefix=/usr/local/share/project/ffmpeg/.pixi/envs/default \
    --enable-gpl --enable-version3 --enable-libdrm --enable-rkmpp --enable-rkrga
make -j $(nproc)
make install
```

## Afterword
After executing all the steps above, the virtual environment will have the ffmpeg specifically optimized for Rockchip chips. The video processing speed is more than 10 times faster than the CPU version~

I'm considering turning the above content into a Pixi environment with tasks, and also learning how to package and upload packages using pixi.
