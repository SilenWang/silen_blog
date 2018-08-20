---
title: Gogs部署
date: 2018-08-18 15:51:59
tags: [Gogs, VPS]
categories: Others
---

买个VPS, 一个月百分之一的流量都用不到, 确实挺浪费的...所以还是多搭点东西利用起来吧, 毕竟一个月$5呢....

<!-- more -->

Gogs是一个Go语言实现的轻量级的git托管服务. 其资源占用之小, 甚至可以直接放在树莓派上.

之前虽然搭建了一次了...但是没有做记录, 这次搭的时候又查了半天....

总的来说, 使用docker部署(简单), 所以实际上配置个文件执行几个就好了, docker的安装不另外说明, archwiki中有了. 主要的安装命令:

```bash
docker pull docker.io/gogs/gogs
```

然后准备一个含有下面内容的文件(`docker-compose.yml`), 执行在有这个文件下的目录执行`docker-compose start`, 就可以开始用了, gogs的访问端口是10080

```txt
version: "2"

services:
  gogs:
    image: docker.io/gogs/gogs

    restart: always
    ports:
      - "10022:22"
      - "10080:3000"
    volumes:
      - /var/gogs:/data

```

然后为了免密码推送, 还需要设置一下ssh, 默认开启gogs之后, 内置的ssh是没有开启的, 需要编辑`/var/gogs/gogs/conf/app.ini`修改到如下样子(这是我的配置文件目录, 根据`volumes`后面实际挂载的位置修改)

```ini
[server]
DOMAIN           = http://1XXX.XXX.XXX.1XXX
HTTP_PORT        = 3000
ROOT_URL         = http://localhost:10080/
DISABLE_SSH      = false
SSH_PORT         = 10022
START_SSH_SERVER = true
OFFLINE_MODE     = false
```

剩下的跟github一样设置就好了~