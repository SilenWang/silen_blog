---
title: Chromebook X11 使用中遇到的两个 Bug 及解决方案
date: 2026-05-30 20:00:00
tags:
  - Chromebook
  - ChromeOS
  - Bug
  - Snapdragon
  - GPU
categories: Others
---

接上篇 {% post_link 从Fydetab-Duo到HP-Chromebook-X11-一个设备折腾党的自我修养 [从 Fydetab Duo 到 HP Chromebook X11] %}，HP Chromebook X11 用了一段时间，体验整体不错，但还是遇到了两个比较烦人的 Bug，记录一下解决方案。

<!-- more -->

## Bug 1：Linux 子系统自启动服务宿主机无法访问

### 问题描述

我在 Chromebook 的 Linux 子系统中部署了 openvscode-server，为了方便，给 openvscode-server 写了一个 systemd service，并通过 `chromeos-autostart` 插件让 Linux 子系统在登录时自动启动。

结果发现：每次开机后，宿主机的浏览器访问 `http://127.0.0.1:3000` 都是连不上的，必须手动进 Linux 子系统执行 `systemctl restart openvscode-server` 才能正常访问。

### 原因分析

其实原因完全不明... 没有看到任何相关的讨论，AI也妹有搜到任何类似问题... 而且原来用 FydeOS 的时候没有出过这问题...

### 解决方案

在 autostart 执行的脚本里，加一句重启命令就行了。比如我在 autostart 脚本中增加：

```bash
sudo systemctl restart openvscode-server
```

改完之后重启测试，宿主机直接访问 `http://127.0.0.1:3000` 就正常了。

## Bug 2：Snapdragon 显卡驱动导致的黑屏与色块问题

### 问题描述

HP Chromebook X11 用的是骁龙 7c Gen 2，这颗 SoC 的 Adreno GPU 在 ChromeOS 上的驱动似乎有不少问题，具体表现有两个：

1. **终端 App 中出现大量黑框色块**：无论是系统自带的 Terminal 还是 crosh，滚动时经常会出现黑色的方块/条纹，严重影响阅读。

2. **打开 opencode 后宿主机动鼠标必黑屏**：这个更严重。在终端中通过命令行打开 opencode，只要回到宿主机桌面移动鼠标，屏幕就会闪黑屏。不是偶发，非常容易出现。

### 原因分析

根据 Gemini 的回答，这两个问题都和 GPU 的 2D canvas 渲染有关。ChromeOS 默认使用 GPU 加速 2D 渲染，但骁龙 7c Gen 2 的 Adreno GPU 驱动在 2D canvas 的实现上可能有兼容性问题：

- 终端 App 大量使用 canvas 来渲染文字和背景，滚动的过程中不断触发重绘，驱动处理不过来就出现了黑块。
- opencode 是一个终端 TUI 应用，其渲染方式触发了某种 GPU 状态泄漏或者死锁，导致整个 ChromeOS 的图形栈受到影响，鼠标移动时的重绘请求直接让显示屏黑了。

### 解决方案

关闭 ChromeOS 的 2D canvas 加速。在 Chrome 中访问：

```
chrome://flags/#disable-accelerated-2d-canvas
```

把这个 flag 设置为 **Disable**

重启之后：
- 终端 App 的黑框色块消失了，文字渲染正常。
- 打开 opencode 后移动鼠标也不会黑屏了。

### 代价

这个特性原本是功过GPU加速2D图形渲染的，关闭之后，相关渲染必然是变卡的...

但说实话，比起黑屏完全没法用，卡一点至少还能忍。希望 ChromeOS 后续能修复这个驱动问题（虽然不是那么现实，似乎7c的设备一直各种问题）。

## 后记

选择商业产品本来是为了更好的可用性和更少的问题... Chromebook 尚且依然有问题... 之前 Fydetab Duo 小问题不断... 也突然觉得能接受了...