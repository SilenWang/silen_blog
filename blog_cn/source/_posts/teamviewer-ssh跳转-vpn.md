---
title: teamviewer+ssh跳转=vpn
categories: Script
date: 2018-12-23 23:00:13
tags: ['VPN', 'SSH跳转', 'TeamViewer', 'VirtualBox', 'Cygwin', 'Windows', 'Linux']
---

<!-- 摘要部分 -->
好久不更了, 今天来记录一下自制vpn的升级版~

<!-- more -->
之前因为一些原因笔记本装了win10, 然后又装回了manjaro. 于是又要重新开始配置一些东西, 其中一项就是回家加班要用的vpn.

由于现在的公司是没有搭vpn的, 所以之前一直用的teamviewer. 但是实际使用下来相当的不方便. 因为我在公司电脑上用的时虚拟机, 经常会出现一些莫名其妙的问题, 再就是网络繁忙的时候会出现明显的卡顿.

之前装win10的那段时间偶然发现teamviewer最新的14版本有一个vpn模式, 试用了一下, 是真的有用的...于是就效仿在上一间公司的自制vpn又如法炮制了一遍:

- 材料准备
    1. virtualbox with micro xp
    2. teamviewer 14 
    3. cygwin with openssh
- tip: 本次试了tiny xp的镜像...发现cygwin在里面没法正常用...没办法还是回百度云下回了原来的备份镜像...

- SSH设置: 配置`~/.ssh/config`, 在面设置好跳板就行了, 其中跳到公司工作电脑的那一级借助teamviewer14的vpn, 地址在连接后可从teamviewer中获取. 地址好像是固定的, 因此直接写在配置文件里面就好了.

```
Host vpn_vm
    User Administrator
    HostName 127.0.0.1
    Port 8022
    IdentityFile ~/.ssh/vpn_vm_key

Host Jump
    HostName 123.456.789.0
    Port 8022
    User silen
    IdentityFile ~/.ssh/jump_key
    ProxyCommand ssh Administrator@vpn_vm -W %h:%p

Host Target
    HostName 123.456.789.1
    Port 22
    User silen
    IdentityFile ~/.ssh/target_key
    ProxyCommand ssh silen@Jump -W %h:%p
```

- 这里做了一个2级的跳转, `vpn_vm`是本地的虚拟机, 里面装好了teamviewer进行vpn连接, openssh则负责使本机能通过该虚拟机进行跳转
- Jump是公司的机器, 上面用vobx装了manjaro并启动了sshd服务
- Target时最终的服务器
- 上面ProxyCommand设置了每1级的跳转命令, 同时设置了密钥以实现无密直接登录
- 登陆时本机直接`ssh Target`就可以了, 相当方便, `scp`, `sftp`也直接用就好
