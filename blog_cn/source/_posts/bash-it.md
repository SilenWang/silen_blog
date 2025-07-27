---
title: bash-it
categories: Script
date: 2018-10-29 21:17:23
tags: ['命令行', 'bash', 'plugin']
---

bash-it是受oh-my-zsh启发而创建的bash插件集合项目. 安装上bash-it不仅可以让bash界面变得更炫酷, 也可以让命令行变得更为好用

<!-- more -->

## 安装

根据bash-it首页的安装教程, 很快就能够安装好

```bash
git clone --depth=1 https://github.com/Bash-it/bash-it.git ~/.bash_it
~/.bash_it/install.sh
```

不过选择的时候要注意, 会让你选择是否保留原来的`.bashrc`以及`.bash_profile`, 当然如果不保留这俩文件也不会被删除, 只是会被备份一下, 如果安装后出了什么问题手动换回来就好了.

## 使用

`bash-it -h`可查看所有可使用的命令, `bash-it show plugins`可以查看所有可用的插件, 目前我只用了`base`和`git`, 其他的日后再慢慢了解好了

![display](https://raw.githubusercontent.com/SilenWang/Gallary/master/bash_it.png)
