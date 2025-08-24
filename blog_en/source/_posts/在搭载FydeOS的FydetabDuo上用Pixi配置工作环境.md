---
title: Random Thoughts on FydeTab Duo
categories: Others
date: 2024-09-24 23:30:45
tags: ['pixi', 'fydeos']
---

<!-- Abstract -->
<!-- more -->

When I initially participated in the crowdfunding on [Indiegogo](https://www.indiegogo.com/projects/fydetab-duo--2), I had already been studying and working on a pure Linux computer for at least six years. My primary systems had changed several times: Ubuntu, Linux Deepin (the ancient version of Deepin OS, which was a downstream of Ubuntu at that time), Linux Mint, Manjaro.

Manjaro seemed to be the most problem-free and rolling updated distribution. Yet, it was inevitably hit with dependency issues after upgrades, preventing important work software from opening and thus delaying work.

Therefore, after visiting the FydeOS website, I had the idea of trying again to use a stable, crash-free "netbook" for work. After all, in terms of computational stability and information security, it is a more ideal solution to perform all bioinformatics analysis tasks on servers or in the cloud. The local machine only needs the most basic tools (SSH) to connect to the server and then it should be able to handle documents, spreadsheets, and presentations at the most basic level.

So, with a little bit of support for the new domestic direction, I joined the crowdfunding for Fydetab Duo, and while waiting, I began to learn about the ChromeOS system, and even later, I tried to work on a Chromebook for a full year.

After a long wait, in May 2024, I finally received the 2024 version of Fydetab Duo, an open-source tablet using an RK3588s processor, which officially used the FydeOS system (though I had already obtained a second-hand engineering version from goofish.com two months ago).

To be honest, even without comparing to full Linux distributions, but just comparing to ChromeOS, there are currently some issues that make this open-source tablet's FydeOS uncomfortable to use.

Some of these issues come from the Linux subsystem. The biggest problems are twofold: first, there is a probability that the subsystem will lose network connectivity after waking up. Once lost, I need to restart the network service within the subsystem, as well as the virtual machine. I haven't found a pattern to this problem, and restarting the virtual machine often takes three to five minutes. This is quite annoying because there may be a computational process running in the subsystem, and during the waiting time, I will engage in compilation and environment preparation for various programs, which often require a network connection. If the network suddenly goes haywire... T T

The other part comes from the hardware. I currently own two versions of FydeTab: a red engineering version and a gray 2024 retail version. The touchpad on the magnetic keyboard cover of both versions is identified as a mouse and cannot use a series of gestures in ChromeOS, which is very unnatural for me after a year of using ChromeOS. At the same time, the sensitivity of the touchpad needs adjustment; it is extremely prone to accidental touches, causing the mouse to move erratically while typing - a problem that was not present on the Chromebook 2017.

Nevertheless, the fact that the company is still in operation and that I have received the device is already quite an achievement. I only hope that these issues can be slowly resolved in the future.
