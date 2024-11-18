---
title: ssh no longer supports rsa by default
categories: Others
date: 2021-10-11 21:06:49
tags: ['ssh', 'openssh', 'manjaro']
---

Today, after upgrading Manjaro, the most basic ssh login suddenly had issues, and the prompts were as follows:

```bash
Unable to negotiate with UNKNOWN port 65535: no matching host key type found. Their offer: ssh-rsa,ssh-dss
lost connection

sign_and_send_pubkey: no mutual signature supported
```

Upon investigation, I found out...well, [OpenSSH deems the ssh-rsa encryption method insecure and has directly disabled this key type for logins by default starting from 8.8](https://www.openssh.com/releasenotes.html)...
<!-- 摘要部分 -->
<!-- more -->

I initially thought that Manjaro lags behind Arch by several versions, so I went to the Archwiki for a solution, but it wasn't mentioned yet... After searching for half an hour, I founded that I could add the following lines in ~/.ssh/config to resolve the issue:

```config
Host *
    PubkeyAcceptedKeyTypes +ssh-rsa
    HostKeyAlgorithms +ssh-rsa
```

The first line means it applies to all hosts, the second line adds `ssh-rsa` back into the allowed use range, and the third line specifies that all hosts should use the `ssh-rsa` algorithm's key.
Through actual testing, both lines are required; without the second line, it prompts that there is no `ssh-rsa` type, and without the third line, it prompts `sign_and_send_pubkey: no mutual signature supported`.

That's all~