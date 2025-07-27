---
title: 钉钉Linux版能正常用了
categories: Daily
date: 2022-02-22 22:22:22
tags: ['钉钉', 'Linux', 'AUR', 'Arch Linux', '办公软件']
---

虽然钉钉某种程度上确实是个压榨人的工具, 但是毕竟工作要用到, 与又被压榨又被恶心... 能去掉一个恶心已经算不错了(说的就是你腾讯)
<!-- 摘要部分 -->
<!-- more -->

感谢[FlyInWind](https://aur.archlinux.org/packages/dingtalk-bin)在AUR上给出的解决方案, 钉钉Linux版本在Arch系下也可以正常使用了, 目前我所有要用的组件都正常(聊天, 工作台的审批, 文档, 项目管理, 在线文档分析), 视频会议本周应该会试用到, 如果也正常...可就齐活了

目前开发者应该是已经把FlyInWind提出的方案直接写到构建文件里了, 使用AUR安装的钉钉直接就把`libgtk-x11-2.0.so`等两个文件删除了, 只要按照dbh625说的在启动脚本加入以下代码即可完美使用.

```bash
export XMODIFIERS="@im=fcitx"
export QT_IM_MODULE="fcitx"
export QT_QPA_PLATFORM=xcb
```

期待所有国产常用软件都能有功能完整的Linux版本的一天
