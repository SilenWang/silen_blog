---
title: 想打包个opencode结果又遇坑了
categories: Others
date: 2026-01-31 15:54:18
tags:
---

上次成功打包了我自己的程序，这想来试个别的，最近 opencode 贼火，我也有在用，刚好 conda 上到目前位置也没有，因此想打包一个。

<!-- more -->

## 打包顺利，但是...

opencode 使用 typescript 编写，同时使用比较新的 bun.js 来进行管理，bun.js 目前在 conda 还没有。另外opencode本身已经为了能在各种环境运行，准备了各种安装包，因此本次编译准备根据 rattler-build 的文档，采用 [repackaging](https://rattler-build.prefix.dev/latest/tutorials/repackaging) 的策略来进行，即不进行编译，而是对编译好的包进行解包然后二次打包，类似之前 AUR 中的各种 deb 转 arch 包的情况。

由于 typescipt 是单文件的程序， 其打包也尤其简单，解压缩出可执行文件放到特定目录就行了。


```toml
package:
  name: opencode
  version: 1.1.35

source:
  - if: target_platform == "linux-64"
    then:
      url: https://github.com/anomalyco/opencode/releases/download/v1.1.35/opencode-linux-x64.tar.gz
      sha256: 451f5a36e2875b5540adf55e8cc9e144902b44959a6f31899fc21876b38b31ae
      file_name: opencode.tar.gz
  - if: target_platform == "linux-aarch64"
    then:
      url: https://github.com/anomalyco/opencode/releases/download/v1.1.35/opencode-linux-arm64.tar.gz
      sha256: e7544ae14afb10e75d28a3623b1fd33d60e17f372106665566ca4e085c2b157b
      file_name: opencode.tar.gz

build:
  number: 0
  script: |
    tar xvf opencode.tar.gz
    mkdir -p $PREFIX/bin
    mv opencode $PREFIX/bin/opencode
    chmod 755 $PREFIX/bin/opencode
    
requirements:
  build:  
    - tar

tests:
  - script:
      - opencode -h

about:
  homepage: https://opencode.ai/
  license: MIT
  summary: 'The open source coding agent.'
  description: |
    The open source coding agent.
  repository: https://github.com/anomalyco/opencode
```

但是问题来了，虽然打包非常顺利，但是与之前无依赖的 Go 程序不同，这次的 conda 包生成后，安装到环境内是不能用的，提示段错误。

## 找问题

于是开始了与 AI 的漫长讨论，最终在 AI 的指导下， 使用 `rattler-build`，解 conda 包，发现二进制被修改过（md5sum不一致）。进而了解到，不论是 `rattller-build` 还是 `conda-build`，都有特别的机制，即对二进制文件进行额外的库路径设置。

使用 `readelf -d opencode` 可以看到，相比获得的二进制文件，conda 包里的 opencode 多了一行 `Library rpath` 内容，就是这个改变导致了程序不能用。

```
Dynamic section at offset 0xd588 contains 32 entries:
  标记        类型                         名称/值
 0x000000000000000f (RPATH)              Library rpath: [$ORIGIN/../lib]
 0x0000000000000001 (NEEDED)             共享库：[libc.so.6]
 0x0000000000000001 (NEEDED)             共享库：[ld-linux-aarch64.so.1]
 0x0000000000000001 (NEEDED)             共享库：[libpthread.so.0]
 0x0000000000000001 (NEEDED)             共享库：[libdl.so.2]
 0x0000000000000001 (NEEDED)             共享库：[libm.so.6]
```

根据文档，这个修改 rpath 的行为是可以关掉的：

```toml
package:
  name: opencode
  version: 1.1.35

source:
  - if: target_platform == "linux-64"
    then:
      url: https://github.com/anomalyco/opencode/releases/download/v1.1.35/opencode-linux-x64.tar.gz
      sha256: 451f5a36e2875b5540adf55e8cc9e144902b44959a6f31899fc21876b38b31ae
      file_name: opencode.tar.gz
  - if: target_platform == "linux-aarch64"
    then:
      url: https://github.com/anomalyco/opencode/releases/download/v1.1.35/opencode-linux-arm64.tar.gz
      sha256: e7544ae14afb10e75d28a3623b1fd33d60e17f372106665566ca4e085c2b157b
      file_name: opencode.tar.gz

build:
  number: 0
  script: |
    tar xvf opencode.tar.gz
    mkdir -p $PREFIX/bin
    mv opencode $PREFIX/bin/opencode
    chmod 755 $PREFIX/bin/opencode
    
  prefix_detection:
    ignore_binary_files: true
  dynamic_linking:
    binary_relocation: false

requirements:
  build:  
    - tar

tests:
  - script:
      - opencode -h

about:
  homepage: https://opencode.ai/
  license: MIT
  summary: 'The open source coding agent.'
  description: |
    The open source coding agent.
  repository: https://github.com/anomalyco/opencode
```

但是诡异的事情来了，在修改设置后，包内的二进制程序和打包前的 md5 一致后，安装 conda 包依然有段错误，检测安装到环境的包，发现这个 rpath 的改动依然存在。

由于 pixi 本身比较新，我自己完全检索不到这个问题的相关讨论，于是我开始了为期一周的，与 AI 的 Battle ... 最后的结果是，pixi global install 的机制如此，没有办法解决。

然而，事实并不如此... 我猛然发现，不论我怎么改 Recipe 文件，似乎打出来的包的名称是固定不变的... conda 和 pixi 都有缓存机制，于是我试着把缓存清除掉... 然后，发现能运行了...

虽然浪费了比较多的时间，好在结果是好的... 目前来说，只依赖 AI，果然还是不完全现实的...