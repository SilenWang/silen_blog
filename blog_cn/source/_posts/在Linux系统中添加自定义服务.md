---
title: 在Linux系统中添加自定义服务
categories: Others
date: 2019-10-05 16:00:31
tags: ['systemd']
---

使用Linux系统的时候经常会需要设置一些开机自启动的东西, 对于使用systemd的发行版, 自行编写服务文件并启用是个非常不错的选择.

<!-- 摘要部分 -->
<!-- more -->

服务文件采用类似ini格式的配置文件编写, 一个基本的服务配置文件包含如下内容:

- `Unit`: 定义的部分, 比如说明什么的
- `Service`: 实际执行执行内容, 包括这个服务是啥类型, 启动时要做什么, 结束时要做什么, 重启要做什么之类的
- `Install`: 基本都是以来关系, 比如服务什么时候启动, 或者在XX之前/之后启动

下面是一个v2ray服务的例子:

```ini
[Unit]
Description=Daemon to start V2ray

[Service]
Type=simple
ExecStart=/usr/bin/v2ray -config /home/silen/script/v2ray/v2ray.json

[Install]
Alias=v2ray
WantedBy=default.target
```

这个服务如果开启了, 则默认自动启动(`WantedBy=default.target`), 启动时执行命令`/usr/bin/v2ray -config /home/silen/script/v2ray/v2ray.json`, 也就是开启v2ray代理. 对于服务的其他动作, 比如重启停止之类没有作定义.

写好上面的文件后就可以把这个文件命名为`v2ray.service`放到指定的位置(`/etc/systemd/system`或`/usr/lib/systemd/system`), 然后以root权限使用`systemctl enable v2ray.service`就可以激活服务(开机后根据配置情况启动), 使用`systemctl start v2ray.service`则可立刻启动服务.

需要注意的是, 调用root权限去启动服务, 那么执行用户也自然是root, 这在某些情况下并不是一个好选择. 因此我之前的同步服务和这次的v2ray服务实际都是以个人用户的权限去开启的. 也就是将上面的服务文件放到`/home/<user>/.config/systemd/user`下, 然后直接`systemctl start --uesr v2ray.service`来启动服务.

这种情况下的服务只会在我当前用户登入后才开始执行, 并且也不需要root权限(用户需要不需要在root group中我倒不清楚).

示例中使用的v2ray程序属于占着前台不放的类型, 因此可以直接使用`simple`类型. 除了该类型还有: `simple`, `forking`, `oneshot`, `notify`, `dbus`, `idle`

当然我自定义的服务一般也就很简单, `sample`和`oneshot`应该够用了.

至于`ExecStart`, `ExecStop`, 和`ExecReload`, 前者是必须的, 后两者则可根据启动命令的情况选择是否需要编写.

最后是`Install`里的`WantedBy=default.target`这个根据我的尝试是必须加的, 因为不加服务不会自动启动...不知道是不是只有user层面的服务才会这样