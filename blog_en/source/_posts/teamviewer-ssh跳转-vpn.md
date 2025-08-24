---
title: teamviewer+ssh jump = VPN
categories: Coding
date: 2018-12-23 23:00:13
tags: ['ssh', 'teamviewer']
---

<!-- Abstract part -->
It's been a while since the last update. Today, I'll record the upgraded version of my homemade VPN!

<!-- more -->
Previously, due to some reasons, my laptop was running Windows 10, and then I switched back to Manjaro. So, I had to reconfigure some things again, one of which was the VPN I used for working from home.

Since our company doesn't have a VPN set up, I've always been using TeamViewer. However, it's quite inconvenient. Because I use a virtual machine on my company computer, there are often mysterious issues, and during busy times, it shows noticeable lag.

During the time when I was running Windows 10, I偶然 discovered that the latest version of TeamViewer (version 14) has a VPN mode, so I tried it out, and it really works... So, I followed the homemade VPN setup from my previous company to create one again:

- Materials needed
    1. VirtualBox with Micro XP
    2. TeamViewer 14 
    3. Cygwin with OpenSSH
- Tip: This time, I tried the Tiny XP image... found that Cygwin couldn't work properly inside it... so I had to download the original backup image from Baidu Cloud...

- SSH Settings: Configure `~/.ssh/config`, set up the jump point there, where jumping to the company's work computer is done via TeamViewer 14's VPN. The address can be obtained after connecting through TeamViewer. It seems to be fixed, so I directly wrote it into the configuration file.

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

- Here, a two-level jump is set up. `vpn_vm` is the local virtual machine where TeamViewer is installed for VPN connection and OpenSSH is responsible for allowing this machine to jump through the virtual machine.
- `Jump` is the company's machine, which has Manjaro installed in VirtualBox and starts the SSHD service.
- `Target` is the final server.
- The `ProxyCommand` sets up the command for each level of jump, and also sets up keys to achieve passwordless login.
- When logging in, you can directly use `ssh Target`, which is very convenient. `scp` and `sftp` can be used directly as well.
```