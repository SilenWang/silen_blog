---
title: Setting up HTTPS using Caddy
categories: Others
date: 2019-05-07 21:59:51
tags: ["caddy", "https"]
---

Previously, I set up a Gogs server but didn't configure HTTPS; it was directly using HTTP. Last time, I tried and found that intercepting passwords was too easy... So I thought about getting a domain name and setting up HTTPS...

<!-- Abstract part -->
<!-- more -->

Firstly, to use HTTPS, you need an SSL certificate. To get a certificate, it seems that free SSL certificates basically require a domain name. So I got the cheapest one from [Namesilo](www.namesilo.com/‎). Domain registration and binding IP can be found in this blog post: [this article](https://zhuanlan.zhihu.com/p/33921436).

After getting the domain, the next step is to set up a reverse proxy server. I chose Caddy because I read that it's simple and easy to use, with automatic application and renewal of SSL certificates. Maybe my luck wasn't good... Although setting it up was very simple, I couldn't get the SSL certificate to work. In the end, I found a project on GitHub called [acme.sh](https://github.com/Neilpang/acme.sh) to apply for SSL certificates. This script is very easy to use, and the author's usage instructions are detailed; following them should be fine. Of course, there's also a more user-friendly tutorial: [this guide](https://www.lizi.tw/soft/5049.html).

Using `acme.sh` will generate four files, among which `fullchain.cer`, `your.domain.key` are the certificate and key.

Since I deployed Gogs and Caddy using Docker, to make things convenient, I wrote both in a single compose file:

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

In the above file, two containers are configured. Gogs only exposes port 22 for SSH service; HTTP ports are not exposed. Caddy can set access directly in its configuration (although I don't know how to do it).
For the Caddy container, note that setting the `ACME_AGREE=true` environment variable is necessary; otherwise, acme will prompt you that you cannot use it... Also, specify the Caddy configuration file `Caddyfile` and the certificates and keys generated above. The port exposure in the Caddy container's configuration seems to be incorrect; I'll look into how to fix it later.

The specific content of the Caddyfile is as follows:

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
```