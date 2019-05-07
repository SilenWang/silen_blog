---
title: 使用Caddy设置https
categories: Others
date: 2019-05-07 21:59:51
tags: ["Caddy", "https"]
---

之前自建了一个gogs, 但是并没有弄https, 直接用的http, 上次自己试了一下发现截密码太容易了...所以思考了一下...还是弄个域名然后上了https...

<!-- 摘要部分 -->
<!-- more -->

首先, 使用https需要申请SSL证书, 然后要申请证书的话, 似乎免费的SSL证书是基本都是要有域名的, 所以我找了最便宜的[Namesilo](www.namesilo.com/‎
)申请了一个. 域名申请和绑定IP见[这个博文](https://zhuanlan.zhihu.com/p/33921436).

有了域名之后下一步就是设置反代理服务器了, 我选择了简单易用的Caddy, 因为我查到的资料都说Caddy设置简单, 易用性完爆Nginx还能自动申请和更新SSL证书. 也许是我脸比较黑...设置虽然很简单完成了, 但是申请SSL证书怎么也用不了. 无奈只有在github上又找了个专门申请SSL证书的项目[acme.sh](https://github.com/Neilpang/acme.sh). 这个脚本使用非常简单, 作者的使用说明也写的详细, 照它说的做就可以. 当然这里也有个[更傻瓜化的教程](https://www.lizi.tw/soft/5049.html)

使用`acme.sh`会生成四个文件, 其中`fullchain.cer`, `your.domain.key`就是证书和密钥.

由于我使用docker部署gogs和caddy, 为了方便把两者都写到一个compose文件里:

```
version: "2"
services:
    gogs:
        image: 08fa8ee955da
        container_name: gogs
        restart: always
        ports:
            - "3022:22"
        volumes:
            - /var/gogs:/data
    caddy:
        container_name: caddy
        image: abiosoft/caddy
        environment:
          - ACME_AGREE=true
        volumes:
          - "~/gogs/Caddyfile:/etc/Caddyfile"
          - "/root/.acme.sh/your.domain/fullchain.cer:/root/fullchain.cer"
          - "/root/.acme.sh/your.domain/your.domain.key:/root/your.domain.key"
        ports:
          - "8080:2015"
          - "80:80"
          - "443:443"
        restart: always
```

上面文件中配置了两个容器, gogs只暴露22端口用于ssh服务, http的端口不暴露, Caddy可以在配置文件里直接设置访问(虽然我不知道怎么做到的).
Caddy的容器则注意设置`ACME_AGREE=true`这个环境变量的设置, 不然acme会提示不能使用... 再就是要指定Caddy的配置文件`Caddyfile`和上面生成的证书和密钥了. Caddy容器的端口暴露这里设置的有问题...之后再看怎么改.

Caddyfile的具体内容如下:

```
your.domain {
    proxy / gogs:3000 {
        header_upstream Host {host}
        header_upstream X-Real-IP {remote}
        header_upstream X-Forwarded-For {remote}
        header_upstream X-Forwarded-Proto {scheme}
    }
    log /var/log/caddy.log
    gzip
    tls /root/fullchain.cer /root/your.domain.key
}
```