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

我在 Chromebook 的 Linux 子系统中部署了 openvscode-server，为了方便，给 openvscode-server 写了一个 systemd service，并通过 `autostart` 插件（`chromeos-autostart`）让 Linux 子系统在登录时自动启动。

结果发现：每次开机后，宿主机的浏览器访问 `http://penguin:3000` 都是连不上的，必须手动进 Linux 子系统执行 `systemctl restart openvscode-server` 才能正常访问。

### 原因分析

这个问题其实跟 openvscode-server 本身无关。Chromebook 的 Linux 子系统（Terminal）在启动时有一个固定的启动顺序：

1. 先启动容器（penguin）
2. 加载用户会话
3. 执行 autostart 插件

但是 autostart 插件触发的时候，systemd 可能还没完全就绪，或者网络相关的服务还没完成初始化。导致 service 虽然被 mark 为 enabled，但实际上并没有成功启动（或者启动时端口还没绑定到容器网络接口上）。

### 解决方案

在 autostart 执行的脚本里，加一句重启命令就行了。比如我在 autostart 脚本中增加：

```bash
systemctl --user restart openvscode-server
```

注意一定要用 `--user`，因为 openvscode-server 的 service 是用户级别的（user service），不是系统级别的。如果不加 `--user`，systemctl 会报错说找不到这个 service。

完整的 autostart 脚本大概长这样：

```bash
#!/bin/bash
# 启动其他服务...
systemctl --user restart openvscode-server
```

改完之后重启测试，宿主机直接访问 `http://penguin:3000` 就正常了。

## Bug 2：Snapdragon 显卡驱动导致的黑屏与色块问题

### 问题描述

HP Chromebook X11 用的是骁龙 7c Gen 2，这颗 SoC 的 Adreno GPU 在 ChromeOS 上的驱动似乎有点小问题，具体表现有两个：

1. **终端 App 中出现大量黑框色块**：无论是系统自带的 Terminal 还是 crosh，滚动时经常会出现黑色的方块/条纹，严重影响阅读。

2. **打开 opencode 后宿主机动鼠标必黑屏**：这个更严重。在终端中通过命令行打开 opencode（一个基于 TUI 的 AI 编程工具），只要回到宿主机桌面移动鼠标，屏幕就会闪黑屏。不是偶发，是必现。

### 原因分析

这两个问题都和 GPU 的 2D canvas 渲染有关。ChromeOS 默认使用 GPU 加速 2D 渲染，但骁龙 7c Gen 2 的 Adreno GPU 驱动在 2D canvas 的实现上可能有兼容性问题：

- 终端 App 大量使用 canvas 来渲染文字和背景，滚动的过程中不断触发重绘，驱动处理不过来就出现了黑块。
- opencode 是一个终端 TUI 应用，其渲染方式触发了某种 GPU 状态泄漏或者死锁，导致整个 ChromeOS 的图形栈受到影响，鼠标移动时的重绘请求直接让显示屏黑了。

### 解决方案

关闭 ChromeOS 的 2D canvas 加速。在 Chrome 中访问：

```
chrome://flags/#disable-accelerated-2d-canvas
```

把这个 flag 设置为 **Enabled**（注意名字是 `disable-accelerated-2d-canvas`，Enabled 表示禁用了加速），然后重启浏览器。

![disable-2d-canvas-flag](https://via.placeholder.com/600x200?text=disable-accelerated-2d-canvas)

重启之后：
- 终端 App 的黑框色块消失了，文字渲染正常。
- 打开 opencode 后移动鼠标也不会黑屏了。

### 代价

不过天下没有免费的午餐，关闭 GPU 加速 2D canvas 的后果就是：**所有浏览器页面明显变卡了**。特别是那些重度使用 canvas 的页面（比如 Jupyter Notebook、VS Code Web、Figma 等），操作起来能感觉到明显的掉帧和延迟。

但说实话，比起黑屏完全没法用，卡一点至少还能忍。希望 ChromeOS 后续能修复这个驱动问题，或者骁龙这边能出一个更新的 GPU 驱动。

## 总结

这两个 Bug 都不算什么大问题，找到了原因之后解决起来也挺快。只是第一次遇到的时候确实让人摸不着头脑，特别是第二个黑屏问题，花了我好几天才定位到是 2D canvas 加速的锅。

Chromebook X11 作为一台二合一设备，硬件底子其实不错，但终究是小众平台，这些边缘 Case 需要用户自己踩坑自己填。记录下来，希望能帮到遇到同样问题的朋友。
