---
title: Mac 与 Windows 文件共享排查记——顺便吐槽一下 Windows 的"屎山"设计
categories: Others
date: 2026-08-10 12:00:00
tags: ['macOS', 'Windows', 'SMB', '文件共享', '网络排查', 'Win11', '吐槽']
---

家里有一台 Mac 一台 Windows，想把文件共享打通。Mac 这边明明已经开了"文件共享"，可 Windows 那边死活看不见机器，也看不见共享目录。折腾了大半天，最后靠着一堆命令行才找到真正的坑。这篇把排查过程记下来，顺便狠狠吐槽一下 Windows 那套从 Win10 到 Win11 越改越乱的"屎山"设计。

<!-- more -->

先说结论：这种问题九成出在两边没对齐的地方——要么是 Windows 的"网络发现"没开，要么是两边的 SMB 协议版本没对上，再要么是凭据没填对。但真正坑人的是，这些设置在 Windows 上散落在好几个地方，找起来能把人逼疯。

## 先确认网络到底通不通：IP 直连

排查的第一步不是看那些花里胡哨的开关，而是先确认底层网络是通的。方法很简单：在 Windows 上 `ipconfig` 查本机 IPv4 地址，然后在 Mac 访达里按 `Command+K`，输入 `smb://IP` 直接连。

```
ipconfig        # Windows 上查 IPv4 地址
# Mac 访达: Cmd+K → smb://192.168.x.x
```

能连上，说明问题基本就出在"网络发现"层面；连不上，那就得往 SMB 协议、凭据这些方向深挖了。

## Windows 侧：网络发现与防火墙

网络发现是 Windows 主动"看见"局域网里其他机器的开关，位置在：

> 控制面板 > 网络和共享中心 > 更改高级共享设置

确认"网络发现"和"文件和打印机共享"都开启了。另外 Windows Defender 防火墙要在"允许应用通过防火墙"里放行"文件和打印机共享"，否则 Mac 那边共享开了也白搭。

## 工作组对齐

Windows 默认的工作组叫 `WORKGROUP`，如果 Mac 那边的工作组跟它对不上，也会互相看不见。去 Mac 上改：

> 系统设置 > 通用 > 共享 > 文件共享 > 详细信息 > WINS 标签页

把工作组名改成和 Windows 一致。

## SMB 协议版本

macOS 10.13 之后默认禁用不安全的 SMBv1，如果对面是 Win7 这种老系统，就得确保启用的是 SMBv2/v3，老协议根本谈不拢。

## 时间同步

两台机器时间差太多，可能会直接导致身份验证失败，看着像"凭据错误"，其实是时钟漂了。顺手把自动同步时间开了。

## 凭据问题：到底该填什么？

Windows 弹窗让你输"网络凭据"的时候，要填的是**被访问机器上的本地账户**，不是微软账户的邮箱，也不是 PIN。常见几种格式：

- `电脑名\用户名`——Windows 上用 `whoami` 查看，比如 `DESKTOP-ABC123\zhangsan`
- `MicrosoftAccount\邮箱`——用微软账户登录时的格式
- `.\用户名`——表示本机用户
- `guest`——无密码共享时用

之前填错过的凭据会一直留着，害人反复失败。去控制面板清一下：

> 控制面板 > 用户账户 > 凭据管理器 > Windows 凭据

删掉和 Mac 相关的记录再重连。

## 关键判断：Mac 能访问 Windows，但 Windows 看不见 Mac

折腾到这一步有个很关键的确认：Mac 能正常访问 Windows 的共享，反过来 Windows 却看不见 Mac。这就基本锁定了——问题出在网络发现，而不是 SMB 本身。

反向访问的方法是直接在 Windows 文件资源管理器地址栏输入 `\\Mac的IP`。注意开头是两个反斜杠，比如 `\\192.168.1.5`。这是绕过网络发现、直接连 IP，只要这一步能通，就说明 Mac 的共享服务本身没问题。

## Windows 连不上 `\\IP` 的深入排查

如果连 `\\IP` 都连不上，那就要按下面这些挨个排查：

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

### 关键诊断结果

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

### 图形界面对应的设置

不想用命令行的话，也有图形界面的对应位置：

- **关闭 SMB 签名要求**：`gpedit.msc` > 计算机配置 > Windows 设置 > 安全设置 > 本地策略 > 安全选项，把"Microsoft 网络客户端：对通信进行数字签名（始终）"设为"已禁用"。
- **启用来宾登录**：计算机配置 > 管理模板 > 网络 > Lanman 工作站，把"启用不安全的来宾登录"设为"已启用"。

不过要注意：`gpedit.msc` 只有 Windows 专业版/企业版/教育版才有，家庭版只能用 PowerShell。而且"不安全的来宾登录"有安全风险，微软自己也建议问题解决后关掉——我是解决完就还原了，毕竟共享这种东西，能用就行，别给自己留后门。

## 顺便吐槽：Windows 的"屎山"设计

这次排查让我真切体会到了 Windows 那套设计的混乱程度，从 Win10 到 Win11 不仅没收敛，反而更碎了：

**第一，shell 有两个。** cmd 和 PowerShell，同一个设置有时候得用 cmd，有时候得用 PowerShell，两套命令语法还不一样。这次连 `Set-SmbClientConfiguration` 这种配置命令都得特定去管理员 PowerShell 里跑，普通用户连门都摸不着。

**第二，设置界面有两套。** "控制面板"和"设置"并存，同一个配置项两个地方都有，但路径、名字都不一样。网络共享相关设置一会儿在控制面板的"网络和共享中心"，一会儿又在新的"设置"应用里，版本一更新位置还来回挪。十年前的老教程能对上控制面板，今年的新教程全在讲"设置"，你根本不知道你该去哪个。

**第三，同一个东西在四个地方之间来回跳。** 光一个 SMB 签名的开关，就能在组策略（gpedit.msc）、控制面板、设置、PowerShell 命令之间反复横跳。排查一个问题，得同时会猜 UI 路径、会背命令行、还得会查策略项——这哪里是给人用的系统，分明是给做过项目交接的老开发留的。

**反观 macOS 侧**，设置高度集中：系统设置 > 通用 > 共享，一条链路全搞定，文件共享、选项、权限都在那一页里，翻两下就齐了。两边对比一下，高下立判。

## 小结

这次的坑总结起来就一句话：**网络是通的，服务是正常的，问题全在 Windows 那套散装协议设置上**。如果你也遇到 Windows 看不见 Mac，先走 IP 直连排除网络问题，再检查网络发现和工作组，最后用 `Test-NetConnection` 和 `lsof -i :445` 这类命令把范围锁死，再对症下药。至于 Windows 的"屎山"……习惯就好，毕竟微软的祖传艺能了。
