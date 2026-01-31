---
title: Trying to Package Opencode and Hitting Another Snag
categories: Others
date: 2026-01-31 15:54:18
tags:
---

Last time I successfully packaged my own program, so I wanted to try something else. Recently opencode has been extremely popular, and I’ve been using it too. Since it’s not yet available on conda, I decided to package it.

<!-- more -->

## Packaging went smoothly, but...

Opencode is written in TypeScript and uses the relatively new bun.js for management, which is not yet available on conda. Additionally, opencode itself already provides various installation packages to run in different environments, so this time I planned to follow the rattler‑build documentation and adopt a [repackaging](https://rattler-build.prefix.dev/latest/tutorials/repackaging) strategy—that is, instead of compiling from source, unpacking the pre‑built binary and repackaging it, similar to how some AUR packages convert .deb packages to Arch packages.

Because the program is a single‑file executable, packaging is especially simple: just extract the executable and place it in the appropriate directory.

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

But here’s the problem: although the packaging itself went smoothly, unlike the previous dependency‑free Go program, the generated conda package couldn’t be used after installation—it reported a segmentation fault.

## Investigating the issue

Thus began a long discussion with AI. Eventually, guided by AI, I used `rattler-build` to unpack the conda package and discovered that the binary had been altered (md5sum didn’t match). Further investigation revealed that both `rattler-build` and `conda-build` have a special mechanism that adds extra library path settings to binary files.

Using `readelf -d opencode` you can see that, compared to the original binary, the opencode inside the conda package had an extra `Library rpath` line, and this change made the program unusable.

```
Dynamic section at offset 0xd588 contains 32 entries:
  Tag        Type                         Name/Value
 0x000000000000000f (RPATH)              Library rpath: [$ORIGIN/../lib]
 0x0000000000000001 (NEEDED)             Shared library: [libc.so.6]
 0x0000000000000001 (NEEDED)             Shared library: [ld-linux-aarch64.so.1]
 0x0000000000000001 (NEEDED)             Shared library: [libpthread.so.0]
 0x0000000000000001 (NEEDED)             Shared library: [libdl.so.2]
 0x0000000000000001 (NEEDED)             Shared library: [libm.so.6]
```

According to the documentation, this rpath modification behavior can be turned off:

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

But here’s the weird part: after making this change, the binary inside the package matched the original md5, yet installing the conda package still resulted in a segmentation fault. Checking the package installed in the environment showed that the rpath modification was still present.

Because pixi itself is relatively new, I couldn’t find any discussion about this issue on my own, so I began a week‑long battle with AI… The final conclusion was that pixi’s global install mechanism works this way and there’s no way around it.

However, that wasn’t actually the case… I suddenly realized that no matter how I changed the Recipe file, the generated package name seemed to stay the same… Both conda and pixi have caching mechanisms, so I tried clearing the cache… And then, it worked…

Although I wasted quite a bit of time, at least the outcome is good… For now, relying solely on AI still isn’t entirely realistic…
