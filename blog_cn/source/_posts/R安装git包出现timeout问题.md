---
title: R安装github上的包出现timeout问题解决
categories: Script
date: 2018-09-05 08:57:17
tags: ['R', 'github', '包安装', 'timeout', 'curl', '手动安装']
---

R安装github上面的包有时会出现连接不到github的问题, 记录一下解决方案

<!-- more -->

本质上是加载`curl`包, 通过curl来下载的样子, 这个可能只使用于linux了.

```r
options(download.file.method = "libcurl")
library('curl')
install_github(PACKAGENAME)
```

- 重新安装的包依然不能返回对象, 尝试手动下载包然后更改源代码, 直接从源代码安装

```r
install.packages(path_to_file, repos = NULL, type="source")
```
