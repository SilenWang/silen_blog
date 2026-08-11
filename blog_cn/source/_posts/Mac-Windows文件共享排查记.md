---
title: Mac 与 Windows 文件共享排查记——顺便吐槽一下 Windows 的"屎山"设计
categories: Others
date: 2026-08-10 12:00:00
tags: ['macOS', 'Windows', 'SMB', '文件共享', '网络排查', 'Win11', '吐槽']
---

我想要将公司 Mac Mini 上的文件，通过 samba 协议共享给其他 Windows 电脑。Mac 这边明明已经开了"文件共享"，可 Windows 那边死活看不见机器，也看不见共享目录。折腾了大半天，最后靠着一堆命令行才找到真正的坑。这篇把排查过程记下来，顺便狠狠吐槽一下 Windows 那套从 Win10 到 Win11 越改越乱的"屎山"设计。

<!-- more -->

先说结论：这种问题九成出在两边没对齐的地方——要么是 Windows 的"网络发现"没开，要么是两边的 SMB 协议版本没对上，再要么是凭据没填对。但真正坑人的是，这些设置在 Windows 上散落在好几个地方，找起来能把人逼疯。

## 先Mac侧确认网络到底通不通：IP 直连

排查的第一步不是看那些花里胡哨的开关，而是先确认底层网络是通的。方法很简单：在 Windows 上 `ipconfig` 查本机 IPv4 地址，然后在 Mac 访达里按 `Command+K`，输入 `smb://IP` 直接连。

```
ipconfig        # Windows 上查 IPv4 地址
# Mac 访达: Cmd+K → smb://192.168.x.x
```

能连上，说明问题基本就出在"网络发现"层面；连不上，那就得往 SMB 协议、凭据这些方向深挖了。

反向访问的方法是直接在 Windows 文件资源管理器地址栏输入 `\\Mac的IP`。注意开头是两个反斜杠，比如 `\\192.168.1.5`。这是绕过网络发现、直接连 IP，只要这一步能通，就说明 Mac 的共享服务本身没问题。

## Windows 连不上 `\\IP` 的排查

当然，如果连 `\\IP` 都连不上，就先要按下面这些挨个排查：

**Mac 侧先检查：**

1. Mac 用户必须设置登录密码——SMB 通常不允许空密码账户。
2. 系统设置 > 通用 > 共享 > 文件共享 > "选项..."，确认用来登录的 Mac 用户名被勾选。
3. 确认共享文件夹对该用户（或 Everyone）的权限是"读与写"或"只读"。
4. 用户名格式用 `Mac电脑名\用户名`，比如 `MacBook-Pro\zhangsan`。

**Windows 侧命令行深挖：**

```powershell
ping MacIP
telnet MacIP 445                       # 或 Test-NetConnection -ComputerName <MacIP> -Port 445
sc.exe query lanmanworkstation
net use \\MacIP /user:用户名 密码      # 返回 1312/1326 等具体错误码
secpol.msc                              # 检查"网络访问：本地账户的共享和安全模型"是否为"经典"
```

`net use` 报的错误码很有用，1312/1326 这类数字直接指明了是哪一层的问题。

**Mac 侧命令行深挖：**

```bash
sudo lsof -i :445   # 确认 smbd 在监听
# 重启 SMB 服务：
sudo launchctl unload -w /System/Library/LaunchDaemons/com.apple.smbd.plist
sudo launchctl load -w /System/Library/LaunchDaemons/com.apple.smbd.plist
# 放宽协议协商版本：
sudo defaults write /Library/Preferences/com.apple.smbd.plist SmbMinVersion "SMB_2_00"
sudo defaults write /Library/Preferences/com.apple.smbd.plist SmbMaxVersion "SMB_3_00"
```

Mac 防火墙里"阻止所有传入连接"要关掉，也可以试试用 `\\电脑名.local` 而不是 IP 访问。

这一步出来两个非常关键的证据，直接锁定了问题范围：

1. `sudo lsof -i :445` 显示 `smbd` 进程在 445 端口监听，而且 IPv4、IPv6 都在监听——**Mac 端 SMB 服务完全正常**。
2. Windows 上 `Test-NetConnection` 的结果：`TcpTestSucceeded : True`（RemoteAddress 172.16.1.25，SourceAddress 172.16.1.21，InterfaceAlias WLAN）——**Windows 到 Mac 的 445 端口完全连通，网络层一点问题都没有**。

网络通、服务正常，那剩下的就只能在 SMB 协议版本协商和身份验证/凭据这两层里找。

## 最终解决方案

最终在 Windows 管理员 PowerShell 里执行了两条命令搞定：

```powershell
Set-SmbClientConfiguration -RequireSecuritySignature $false -Force
Set-SmbClientConfiguration -EnableInsecureGuestLogons $true -Force
```

执行完再连，一下就通了。这两条分别是关闭 SMB 签名要求和允许不安全的来宾登录——也就是说，真正卡住问题的其实是在协议协商这一层，而不是什么网络、防火墙。

## 吐槽：Windows 的"屎山"设计

这次排查让我真切体会到了 Windows 那套设计的混乱程度，从 Win10 到 Win11 不仅没收敛，反而更碎了：

**第一，shell 有两个。** cmd 和 PowerShell，同一个设置有时候得用 cmd，有时候得用 PowerShell，两套命令语法还不一样。这次连 `Set-SmbClientConfiguration` 这种配置命令都得特定去管理员 PowerShell 里跑，普通用户连门都摸不着。

**第二，设置界面有两套。** "控制面板"和"设置"并存，同一个配置项两个地方都有，但路径、名字都不一样。网络共享相关设置一会儿在控制面板的"网络和共享中心"，一会儿又在新的"设置"应用里，版本一更新位置还来回挪。十年前的老教程能对上控制面板，今年的新教程全在讲"设置"，你根本不知道你该去哪个。

**第三，同一个东西在四个地方之间来回跳。** 光一个 SMB 签名的开关，控制面板、设置、PowerShell 命令之间反复横跳。排查一个问题，得同时会猜 UI 路径、会背命令行、还得会查策略项——这哪里是给人用的系统，分明是给做过项目交接的老开发留的。

**反观 macOS 侧**，设置高度集中：系统设置 > 通用 > 共享，一条链路全搞定，文件共享、选项、权限都在那一页里，翻两下就齐了。两边对比一下，高下立判。

## 小结

这次的坑总结起来就一句话：**网络是通的，服务是正常的，问题全在 Windows 那套散装协议设置上**。如果你也遇到 Windows 看不见 Mac，先走 IP 直连排除网络问题，再检查网络发现和工作组，最后用 `Test-NetConnection` 和 `lsof -i :445` 这类命令把范围锁死，再对症下药。至于 Windows 的"屎山"……只能自适应了……
