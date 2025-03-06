---
title: clash_nps self-built inner net penetration using VPN
categories: Script
date: 2023-01-08 22:03:20
tags: ['nps', 'clash']
---

Since the start of the pandemic, I have had some remote work needs. Today, I finally put together a relatively convenient and nearly VPN-like experience for accessing inner networks.
<!-- Abstract part -->
<!-- more -->

## Past Used Schemes

My initial requirement was simply to log in to the target server through SSH. Therefore, I initially bought a VPS as a jumpboard machine. On the target machine, I used SSH or autossh to start reverse proxy services, allowing me to access the 22 port of the target machine via the port on the VPS. The biggest advantage of this scheme is its simplicity; SSH is something that all Linux systems come with by default. Additionally, once connected via SSH, you can use it for reverse proxying, theoretically enabling access to any inner network IP and port. However, its drawbacks are that as more access needs arise, each new IP/port requires logging into the target machine again and running a reverse proxy command. This makes it no longer simple or user-friendly... Furthermore, if power is lost and needs to be restarted... All these steps have to be repeated... Extremely troublesome.

Later, IT started enabling NPS (Network Penetration Server), which is an inner net penetration tool written in Go. Installation and usage are extremely simple. After setting up the server and client, if there are new access requirements, you can create a private link on the NPS management page to set it up. If encountering restarts, restarting the NPS client is sufficient. Additionally, with administrative permissions, the NPS client can also be configured to start automatically. This saves a lot of work. However, when too many things need to be accessed, each port still requires running an NPC process locally. Furthermore, because it maps ports to the local machine, inner network addresses become 127.0.0.1, which can cause some services to have issues (e.g., Gitea's content).

Recently, while checking the NPS management interface, I learned that the Socks proxy mode can be used to access various inner net resources and implement VPN functionality. Therefore, I tried it out again and finally achieved accessing inner network resources using the Socks mode and Clash client, without affecting external internet traffic.

## Current Linux Scheme
### Required Software/Hardware

- A VPS with a public IP, referred to as `VPS`, which also needs to be able to configure port forwarding permissions on this VPS.
- A host within the target inner net that can access the `VPS`. This machine is referred to as `jumpboard`, and the request for accessing the inner net actually originates from this machine.
- NPS server/client software, which can be obtained from GitHub [here](https://github.com/ehang-io/nps).
- A Clash client software that can intercept global traffic. For Linux, you can use the Linux version of [Clash For Windows](https://github.com/Fndroid/clash_for_windows_pkg) or [Clash Verge](https://github.com/zzzgydi/clash-verge). Note that there was a permission issue with Clash Verge's Tun mode before, but I'm not sure if it has been fixed.

### Implementation Steps

1. Set up NPS on the `VPS` and start the client on the `jumpboard`. Then configure the Socks proxy settings in the NPS service interface. For specific operations, refer to [this blog](https://blog.csdn.net/ha0shenqi/article/details/111194246). Note down the account and password for the Socks service.
2. Install Clash For Windows and follow the [official documentation](https://docs.cfw.lbyczf.com/contents/tun.html#linux) to enable Tun mode (global mode, which will create a virtual network card to intercept all traffic on your machine).
3. Since I already have a converted configuration file for Clash, open the existing configuration file and manually add the Socks5 service from the `VPS`. Also, fill in the Clash rules yourself. An example configuration segment is as follows (the rules should be written at the beginning because once a match is found, it won't continue matching):

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

- Tip: Additionally, if the rule matches a single address, it should be written as `XXX.XXX.XXX.XXX/32`. If it matches an entire subnet, it should be written as `XXX.XXX.XXX.0/16`. Of course, writing a specific address with `/16` will cause the entire subnet to be proxied.

4. Restart the Clash service to access.

## Current ChromeOS / ChromiumOS / Android Scheme

For these systems, only the Clash client changes. Use the `Clash for Android` app from the store. For ChromeOS, you need to first enable the Android subsystem and then install it from the store (before installation, set the system proxy; after installation, close the system proxy and use Clash client to route traffic). Since [Android's VPN inherently intercepts all traffic](https://github.com/Kr328/ClashForAndroid/issues/1647), you don't need to configure mode settings. Just configure the rules directly. However, by default, the client does not proxy requests to inner network IPs. Find the `Settings > Network > Bypass Private Networks` option and disable it.

## Future Plans

Since using Clash for rule-based traffic diversion, this scheme can actually set up simultaneous access to multiple inner net resources, as long as the IP addresses do not conflict. In the future, I will change my home network's IP segment directly...
