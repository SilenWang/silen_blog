---
title: ssh默认不支持rsa了
categories: Others
date: 2021-10-11 21:06:49
tags: ['ssh', 'openssh', 'manjaro']
---

今天升级manjaro后最基本的ssh登陆突然出问题了, 提示有几种:
```bash
Unable to negotiate with UNKNOWN port 65535: no matching host key type found. Their offer: ssh-rsa,ssh-dss
lost connection

sign_and_send_pubkey: no mutual signature supported
```

一查发现...好嘛, [openssh觉得ssh-rsa加密方式不安全, 直接从8.8开始默认不允许这种密钥用于登陆了](https://www.openssh.com/releasenotes.html)...
<!-- 摘要部分 -->
<!-- more -->

本来想这manjaro落后arch一定版本, 所以又跑去archwiki找解决办法, 没想到还没写到里面去...兜兜转转查了半小时, 发现可以在`~/.ssh/config`里面加这么一段解决:
```config
Host *
    PubkeyAcceptedKeyTypes +ssh-rsa
    HostKeyAlgorithms +ssh-rsa
```

第一行说明对所有主机生效, 第二行是将`ssh-rsa`加会允许使用的范围, 第三行是指定所有主机使用的都是`ssh-rsa`算法的key.
实测两行都得要写才行, 没有第二行提示没有`ssh-rsa`这么个类型, 没有第三行就提示`sign_and_send_pubkey: no mutual signature supported`.

以上~
