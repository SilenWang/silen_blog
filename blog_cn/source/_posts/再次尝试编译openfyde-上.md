---
title: 再次尝试编译openfyde（上）
categories: Script
date: 2025-06-15 10:50:34
tags: ['openfyde']
---

今年年初的时候，就尝试编译openfyde，想着学习下编译和修改openfyde系统，也利用一下我多余的一台工程机。然而半年都过去了... 上次编译出镜像的问题都还没解决...

<!-- more -->

看着最近最新的r132代码也在逐步更新了，那我先... 尝试再编译一次吧... 

本次基本遵循{% post_link 编译openfyde [之前的步骤] %}，即装系统依赖，下载`depot_tools`，克隆代码，进入chroot环境然后开始编译。不同的是，目前 r132 版本的代码应该是还没有完整上传，因此编译中还会碰到一些问题：

- `chromeos-kernel-6_1-6.1.115-r1676`的克隆地址还没有修改，无法获取源代码因此无法编译
- `libv4lplugins-0.0.1`代码还有Bug，需要再加一个补丁修复：
```
Index: libv4l-rkmpp/src/libv4l-rkmpp.c
===================================================================
--- libv4l-rkmpp.orig/src/libv4l-rkmpp.c
+++ libv4l-rkmpp/src/libv4l-rkmpp.c

@@ -1206,5 +1206,4 @@
 	.init = &plugin_init,
 	.close = &plugin_close,
 	.ioctl = &plugin_ioctl,
-	.mmap = &plugin_mmap,
 };
```

目前我的破电脑还在编译着编译着就死机的状态... 回头完成了再更新... 