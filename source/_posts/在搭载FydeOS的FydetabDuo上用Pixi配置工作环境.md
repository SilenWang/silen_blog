---
title: 关于FydeTab Duo的碎碎念
categories: Bioinfomatic
date: 2024-09-24 23:30:45
tags: ['pixi', 'fydeos']
---


<!-- 摘要部分 -->
<!-- more -->

我当初参加[Indintogo上的众筹](https://www.indiegogo.com/projects/fydetab-duo--2)时，已经在纯Linux的电脑上学习/工作至少六年了。主役系统换了又换：Ubuntu、Linux Deepin（现在的Deepin OS的远古版本，当时是Ubuntu的下游）、Linux Mint、Manjaro。

Manjaro一轮下来最不容易出问题，且持续更新的系统了。却也免不了在升级后碰上依赖爆炸，导致重要工作软件打不开，耽误工作。

故，在看过FydeOS的网站后，有了再次尝试用一台稳定、不崩溃的"上网本"来工作的想法。毕竟不论从计算稳定，还是从信息保密角度，将全部的生物信息分析工作放在服务器或者云端，是更理想的解决方案。本地机器只要具备连通到服务器的最基本工具（SSH），然后能最低限度的处理文档、表格和幻灯，便足矣。

于是当时抱着一点点支持新国产方向的想法，参加了Fydetab的众筹。并在等待中开始了解ChromeOS系统，甚至后来真的用Chromebook 2017工作了整整一年。

经历漫长的等待后，在2024年的五月份，我终于拿到2024版本的Fydetab Duo，一台使用RK3588s处理器的开源平板，正式用上了FydeOS系统（其实2个月前已经从咸鱼淘到了二手工程版就是）。

老实来讲，即使不与完整的Linux发行版做比，只单单比较ChromeOS，这块开源平板上的FydeOS目前也存在一些令人难受的问题。

其中一部分来向自于Linux子系统。目前最大的问题有二，一是子系统有概率会在黑屏重点亮后，失去网络连接，一旦失去，我需要重启子系统内的网络服务，以及重启虚拟机。我没找到这个问题的规律，并且重启虚拟机很多时候会需要三五分钟。这就比较恼人，因为子系统内可能正跑着一个计算进程，然后在等待的时候，我会进行一些程序编译和其他环境准备的工作，而这些工作很多时候都需要网络。如果这时候，网络突然抽风了.....T T

另外一部分，就来自于硬件了，我目前拥有两个版本的FydeTab，一个是红色的工程版，一个是灰色的2024正式版。两个版本磁吸键盘套上的触控板都会被识别成鼠标，无法使用ChromeOS下的一系列手势，让已经用了一年ChromeOS的我非常不习惯。同时触控板的灵敏度有待调整，当前及其容易误触，打字的时候鼠标老乱动，同样的问题在Chromebook 2017在不曾出现过。

不过即使如此，公司还在，设备也拿到手了已实属不易，只希望后续能慢慢把这些问题都解决掉了。