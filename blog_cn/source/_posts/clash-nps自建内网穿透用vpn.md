---
title: clash_nps自建内网穿透用vpn
categories: Others
date: 2023-01-08 22:03:20
tags: ['内网穿透', 'VPN', 'socks5代理', '远程办公', '流量代理', 'nps', 'clash', 'socks5']
---

自疫情开始, 我就多少有一点远程办公的需求, 今天终于组合处了一个相对方便, 近似内网vpn使用体验的方案.
<!-- 摘要部分 -->
<!-- more -->

## 过去使用过的方案

我最开始的需求很简单, 只是通过ssh登陆目标服务器就好, 因此最开始的直接买了个VPS, 用这个VPS做跳板机, 在目标机器上用ssh或者autossh开启反向代理, 使得我可以经由VPS下的端口访问目标机器的22端口, 从而远程登陆目标机器. 这个方案最大的好处就是简单, ssh是所有linux系统必然会带的东西, 同时ssh能连通后, 可以通过ssh进行反代理, 因此理论上也能实现访问任意的内网ip和端口. 不过其弊端是, 当需要访问的内容越来越多, 每访问一个ip/端口就要通过跳板登陆目标机器, 在该机器上再运行一条反代理命令, 就与简单易用没什么关系了... 而且万一断电要重启... 这些东西都要再来一次... 极其的麻烦

之后IT的同时开始启用NPS, 这是个用Go编写的内网穿透工具, 安装使用极其简单, 在设置好服务端和客户端后, 如果有新的访问需求, 在NPS自带的管理页面作相应设置新建一个私密链接就好, 如果遇到重启, 只要重启NPS的客户端就好, 另外如果有管理员权限, NPS客户端还能自己设置自启动服务. 这样就节省了很大一部分工作, 不过当要访问的东西实在是太多, 每一个端口还是要在本地运行一个NPC进程, 且因为是端口映射到本地, 内网地址会变成127.0.0.1, 有些服务访问起来会出问题(Gitea的部分内容就是)

前段时间在查看NPS管理界面时, 界面介绍Socks代理模式可用于访问各种内网资源, 实现VPN的功能, 于是我就又尝试了一下, 最后靠Socks模式和Clash客户端实现了内网资源的访问, 同时又不影响访问外网的流量

## 当前Linux下的方案
### 需要准备的软/硬件

- 一台有公网IP的VPS, 后面简称`VPS`, 同时需能够配置这台VPS的端口开放权限
- 一台位于目标内网中的, 能够访问到`VPS`的主机(系统没太大限制, 因为NPS的客户端是跨平台的), 改机器后面简称`跳板`, 访问内网的请求实际上是到达`跳板`后发出的
- NPS服务端/客户端软件, 可从github[获取](https://github.com/ehang-io/nps)
- 可以进行截取全局流量的Clash客户端软件, Linux下可以使用[Clash For Windows](https://github.com/Fndroid/clash_for_windows_pkg)的Linux版本或者[Clash Verge](https://github.com/zzzgydi/clash-verge), 不过Clash Verge之前开Tun模式似乎存在权限问题, 不知道目前修复了没有

### 实施方法

1. 在`VPS`搭建好NPS, 在`跳板`上开启VPS的客户端, 然后在NPS服务界面进行Socks代理设置, 具体操作方式可参见[这个博客](https://blog.csdn.net/ha0shenqi/article/details/111194246), 注意记录下socks服务的账户密码
2. 安装好Clash For Windows, 然后按照[官方文档](https://docs.cfw.lbyczf.com/contents/tun.html#linux)开启Tun模式(全局模式, 将创建虚拟网卡截取本机的所有流量)
3. 由于我已经通过转换服务产生了Clash配置文件, 因此需要打开已有的配置文件, 手动添加`VPS`上的Socks5服务, 同时自行填写访问内网资源的Clash规则, 示例配置段如下(规则中的rules要写在最前, 因为前面匹配到了后面就不继续匹配了):

```yaml
proxies:
  - name: NPS_SOCKS
    server: 11.11.11.111
    port: 2233
    type: socks5
    username: socks5
    password: socks123456


proxy-groups:
  - name: SOCKS5
    type: select
    proxies:
      - NPS_SOCKS
      - DIRECT


rules:
  - IP-CIDR,192.168.0.111/32,NPS_SOCKS
```

- tips: 另外, ip规则如果匹配单个地址, 写法为`XXX.XXX.XXX.XXX/32`, 如果匹配整个网段, 则为`XXX.XXX.XXX.0/16`, 当然如果写了个具体地址配上`/16`, 整个网段都会被代理

4. 重启`Clash`服务即可访问


## 当前Chromeos / ChromiumOS / 安卓下的方案

这些系统下只有使用的Clash客户端发生变化, 用商店里的`Clash for Android`就好, chromeos需要先开启安卓子系统然后从商店安装(安装前设置系统代理, 之后关闭系统代理通过Clash客户端走流量), 由于[安卓的VPN本来就是截取所有流量的](https://github.com/Kr328/ClashForAndroid/issues/1647), 所以不需要作模式的设置, 直接配置好规则就行, 不过客户端默认设置下内网ip段的请求是不走代理的, 找到`设置 > 网络 > 绕过私有网络`选项, 把选项关闭即可

## 后续

由于使用了Clash进行规则分流, 所以这套方案其实可以设置同时访问多个内网的资源, 当然前提是ip不能冲突, 之后把家里的ip段直接改掉好了...
