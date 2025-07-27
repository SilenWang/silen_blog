---
title: 使用Syncthing进行文件同步
categories: Others
date: 2019-08-18 00:03:33
tags: ['文件同步', '跨设备同步', 'Syncthing', 'systemd', '开源工具']
---

只要使用多台设备进行工作, 就一定会涉及到不同设备之前互相同步的问题. 现在已经不是七八年前云服务刚兴起, 便宜又大碗的时代了. 一个实惠好用又不会突然跑路或改用户协议的服务商怕是并不存在, 所以还是要靠自己来搭了...
<!-- 摘要部分 -->
<!-- more -->

Syncthing在[小众软件](https://www.appinn.com/syncthing/)和[异次元软件](https://www.iplaysoft.com/syncthing.html)都有介绍, 是一款使用简单且开源的同步软件.

由于官方已经提供了打包好的可执行文件, 三大平台都是下载运行文件后即可使用, 因此我就不再重复了. 这里主要记录如何在我的电脑和VPS上将其设置为服务来运行.

CentOS7和Manjaro都有systemd, 且操作命令都是`systemctl`, 所以操作一样的. 

首先找到Syncthing执行文件包, 里面的`Syncthing/etc/`下有已经准备好的各种服务配置文件实例. 我使用的是`Syncthing/etc/linux-systemd/user/syncthing.service`. 打开这个文件, 将其中的`[Service]`部分的`ExecStart`项目中的程序路径改成自己的程序位置. 我是将这个服务作为当前用户的服务而不是系统级服务使用, 所以将这个文件复制到`~/.config/systemd/user/`, 然后运行命令启动服务即可:

```bash
systemctl --user enable syncthing.service
systemctl --user start syncthing.service
```

启动后可像其他服务一样用`ststemctl`来查看服务运行状态.

![Syncthing](https://raw.githubusercontent.com/SilenWang/Gallary/master/Syncthing.png)
