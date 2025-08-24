---
title: R Install GitHub Package Timeout Issue Solution
categories: Coding
date: 2018-09-05 08:57:17
tags: ['R']
---

When installing packages from GitHub in R, sometimes you may encounter a connection timeout issue with GitHub. Here's how to solve it.

<!-- more -->

Essentially, this involves loading the `curl` package and using `curl` to download the package. This might only work on Linux systems.

```r
options(download.file.method = "libcurl")
library('curl')
install_github(PACKAGENAME)
```

- If the installed package still cannot return an object, try manually downloading the package and then modifying the source code to install directly from the source.

```r
install.packages(path_to_file, repos = NULL, type="source")
```
```