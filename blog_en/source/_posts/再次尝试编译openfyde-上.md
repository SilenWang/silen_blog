---
title: Trying to Compile openfyde Again (Part 1)
categories: Script
date: 2025-06-15 10:50:34
tags: ['openfyde']
---

At the beginning of this year, I tried to compile openfyde, hoping to learn about compiling and modifying the openfyde system, and also to make use of my extra engineering machine. However, half a year has passed... and the issues from the last compilation still haven't been resolved...

<!-- more -->

Seeing that the latest r132 code is being gradually updated recently, I decided to... try compiling it once again...

This time I basically followed the [previous steps](编译openfyde.md): installing system dependencies, downloading `depot_tools`, cloning the code, entering the chroot environment and starting compilation. The difference is that the r132 version code hasn't been fully uploaded yet, so there are still some issues during compilation:

- The clone address for `chromeos-kernel-6_1-6.1.115-r1676` hasn't been modified yet, making it impossible to get the source code and compile
- There are still bugs in the `libv4lplugins-0.0.1` code that require an additional patch to fix:
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

Currently my old computer keeps crashing during compilation... I'll update this post when it's done...
