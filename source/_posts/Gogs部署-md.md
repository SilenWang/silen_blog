---
title: Gogs部署
date: 2018-09-21 18:32:00
tags: [Gogs, VPS]
categories: Others
---

买个VPS, 一个月百分之一的流量都用不到, 确实挺浪费的...所以还是多搭点东西利用起来吧, 毕竟一个月$5呢....

<!-- more -->

# 部署步骤
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
DOMAIN           = http://XXX.XXX.XXX.XXX
HTTP_PORT        = 3000
ROOT_URL         = http://localhost:10080/
DISABLE_SSH      = false
SSH_PORT         = 10022
START_SSH_SERVER = true
OFFLINE_MODE     = false
```

剩下的跟github一样设置就好了~

# 无法添加多个密钥问题解决(update@20180921)

gogs 虽然自带有图形化的密钥设置界面, 但是可能时我当时拉取的版本有bug, 添加第一个密钥时是没问题的, 但是当有多台机器, 需要添加一个以上的密钥时, 在图形化界面虽然显示为添加成功, 但并不能正常的通过ssh进行免密登陆, 还需要进一步设置:

- 因为我的docker容器是把`/data`挂载到容器外目录的, 所以不需要登陆docker内系统直接进行修改.
- 对应的文件是:`挂载目录/gogs/git/.ssh/authorized_keys`. 自行对照格式添加公钥内容进去就可以, 比如:
```txt
command="/app/gogs/gogs serv key-1 --config='/data/gogs/conf/app.ini'",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-rsa AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA test
```

- 保存后即可登陆, 不需要重启容器