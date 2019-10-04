---
title: 台式机置备计划(下)
categories: Daily
date: 2019-09-08 12:35:13
tags: ['GPU Passthrough', 'Majaro', 'kvm']
---

在千辛万苦搞定了硬件之后, 接下来就是软件了.

<!-- 摘要部分 -->
<!-- more -->

按照我之前的规划, 我的机器需要达到一机多用的效果: 宿主机linux用来做计算, 虚拟机windows用来打游戏.

但凡用过虚拟机的人应该都知道...虚拟机的性能一直是个的问题, 就算牛X如vmware和mac专属的parallel, 其图形性能也还是不如原生硬件. 因此也就有了PCI Passthrough这种技术, 即将PCI插槽上的设备从宿主机屏蔽, 直接分配给虚拟机使用. 获得了实际的硬件, 没有中间商, 性能也就几乎与原生无异了. 这也是我当初买了个有两条显卡PCI插槽主板的原因, 因为如果日后要虚拟出两台打游戏用的电脑, 必须给每台机器都分配一张显卡, 然后宿主机启动也是要显卡的, 所以CPU需要有内建显卡.

本次虚拟化配置我使用的宿主系统为Manjaro, 参考的文档/视频包括:

- [ArchWiki上的GPU Passthrough教程](https://wiki.archlinux.org/index.php/PCI_passthrough_via_OVMF_(%E7%AE%80%E4%BD%93%E4%B8%AD%E6%96%87))
- [Linus Tech Tips的虚拟MacOS视频](https://www.bilibili.com/video/av54526748)
- [一个关于promox的错误讨论帖]()

由于这份配置是在基于Arch的发行版上进行的, 因此具体设置用的命令和组件会和中文网站中搜索到的材料比较不一样(基本都是基于Ubuntu的).

## 操作步骤

### 1. 准备系统盘

- 宿主系统: manjaro-architect-18.0.4-stable, 基于Arch的衍生发行版, 特点是简单易上手, 虽然被Arch教派认为是邪道, 单阻止不了其越来越高的人气(甚至我昨天看到Manjaro的社区和Ubuntu一样开了个公司). 使用architect而不是预搭载桌面的livecd主要是因为镜像小, 所有软件包全程联网下载就好.
- 子系统: 与宿主系统不一样, 子系统除了准备一份要安装到宿主机中的windows之外, 可能还需要准备一份**支持EFI启动的ISO镜像**. 因为ArchWiki上的教程使用了OVMF(这会让直通更简单, 少一些在虚拟BIOS上配置的过程), 用了这东西之后虚拟机内使用的就是EFI(only), 而不是BIOS了. 这就要求安装系统的ISO镜像一定要支持EFI启动. 但是我下载了很多个系统安装镜像都不支持...因此只好专门找了个支持UEFI的PE镜像, 用这个镜像引导启动, 然后安装另外一份安装镜像上的系统.

### 2. 安装宿主系统

使用宿主系统盘在电脑上安装好Manjaro, 详细步骤其他地方找得到, 我这里就先略了, 以后有空再补截图(咕咕咕).

### 3. 进行GPU隔离设置

进行直通就是让虚拟机直接使用被直通的硬件, 因此需要防止宿主机对硬件进行访问. 具体设置为:

- 启动IOMMU
    - 打开配置文件(`/etc/default/grub`)更改内核参数, 以启用IOMMU:
        ```
        ...
        GRUB_CMDLINE_LINUX_DEFAULT="quiet udev.log_priority=3 audit=0 amd_iommu=on iommu=pt"
        GRUB_CMDLINE_LINUX=""
        ...
        ```
    - 进行上面的编辑后需要重新生成生成`grub.cfg:
        ```
        sudo grub-mkconfig -o /boot/grub/grub.cfg
        ```
- 编辑配置文件使用vcfio屏蔽准备隔离的显卡组件
    - 使用`lspci -nnk`命令查看所有的PCI硬件信息, 从中找到需要直通的硬件, 并记录其ID:
        ```
        01:00.0 VGA compatible controller [0300]: Advanced Micro Devices, Inc. [AMD/ATI] Baffin [Radeon RX 460/560D / Pro 450/455/460/555/555X/560/560X] [1002:67ef] (rev e5)
            Subsystem: Tul Corporation / PowerColor Baffin [Radeon RX 460/560D / Pro 450/455/460/555/555X/560/560X] [148c:2385]
            Kernel driver in use: vfio-pci
            Kernel modules: amdgpu
        01:00.1 Audio device [0403]: Advanced Micro Devices, Inc. [AMD/ATI] Baffin HDMI/DP Audio [Radeon RX 550 640SP / RX 560/560X] [1002:aae0]
            Subsystem: Tul Corporation / PowerColor Baffin HDMI/DP Audio [Radeon RX 550 640SP / RX 560/560X] [148c:aae0]
            Kernel driver in use: vfio-pci
            Kernel modules: snd_hda_intel
        ```
    - 上面的信息是我屏蔽之后的信息, 因此两个设备的`Kernel driver in use`都显示为`vfio-pci`, 代表屏蔽成功, 屏蔽前应该是安装系统时安装的开源驱动; 需要记录的ID跟在设备名称的中括号中, 在上面的例子中分别是显示设备`[148c:2385]`和音频设备`[1002:aae0]`(现在的显卡很多都带音频输出...所以会显示为两个设备实际是要一起屏蔽的)
    - 创建`/etc/modprobe.d/vfio.conf`, 并在其中加入`options vfio-pci ids=<id to block>`
    - 编辑`/etc/mkinitcpio.conf`, 在其中找到下面项目并添加如相应内容:
        ```
        ...
        MODULES=(... vfio_pci vfio vfio_iommu_type1 vfio_virqfd ...)
        ...
        HOOKS=(... modconf ...)
        ...
        ```
    - 编辑完成之后, 更改是不会生效的, 需要重新生成`initramfs`:
        - **注意:** 由于我对这个东西是什么并没有概念, 并且ArchWiki上相关的内容我也看不懂, 所以这里记录的是我进行的实际操作, 并不保证是最正确的操作
        - 执行命令重新生成`/boot`下的启动用`img`文件: `sudo mkinitcpio -g /boot/initramfs-4.19-x86_64.img`(建议把原来的备份)
        - 覆盖原来的`img`后, 需要重新启动宿主机, 如果设置成功, 则对应设备的驱动应该为`vfio-pci`
        - 在进行覆盖后, 宿主系统无法访问被屏蔽的硬件, 如果有需要访问的, 可以在启动时在grub菜单选择通过`initramfs-4.19-x86_64-fallback.img`启动, 这个没有被改动过, 可以正常访问硬件.

### 4. 安装并配置虚拟机

在做好硬件的设置后, 下面就要进行软件的配置了.
- 安装虚拟机管理软件
    + tip: 这里安装的是使用虚拟机的软件QEMU和简化QEMU设置的VirtManager和libvirt, kvm本身是linux内核中的一个模块, 只要在主板中设置开启虚拟化就可以使用了, 无须另外安装什么
    ```bash
    sudo pacman -S qemu libvirt virt-manager ovmf
    ```
- 创建虚拟机
    - 打开virt-manager, 点击创建虚拟机
        ![add_vm.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/add_vm.jpg)
    - 选择安装镜像(非启动用), 由于下载的镜像不是官方的, 无法自动识别是啥系统, 所以去掉自动检测然后手动选择`Win10`然后进入下一步
        ![select_iso.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/select_iso.jpg)
    - CPU和内存分配
        ![cpu_mem.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/cpu_mem.jpg)
    - 选择创建镜像(没啥特殊图略...)
    - 确认信息, 选择网络为桥接, 并选择进一步编辑
    - 分配需要直通的设备
        ![pci_device.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/pci_device.jpg)
    - 由于我手头的安装镜像不支持UEFI启动, 因此另外找了一个支持但是里面没有安装文件的启动镜像, 需要对启动顺序做点小设置, 保证挂载启动镜像并从它启动
        ![boot.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/boot.jpg)
    - 编辑CPU信息, 主要是勾选`复制主机CPU配置`(等同于`host-model`)
        ![set_cpu.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/set_cpu.jpg)
    - 编辑使用的BIOS为OVMF(具体参照ArchWiki, 图忘截了...)
    - 点击左上角的开始安装正常安装系统即可

### 5. 音频设备BUG修复及其他

在完成上述步骤安装好系统后, 一台有物理显卡的虚拟机就搞定了, 将显卡的显示接口接到任意显示设备就能快乐的玩游戏了. 但是这时依然存在一个问题, 虽然显卡是自带音频设备的, 理论上一条HDMI线就可以把音频传到我显示器的扬声器上播放才是. 但是实际上在虚拟机里面连没有找到显卡附带的音频设备. 这个问题苦恼了我一个星期, 最后在Promox(一个专门用来作虚拟化宿主的发行版)的论坛上看到了类似的问题, 说这个问题是Q35 4.0的一个Bug导致的, 可以在虚拟机配置文件中将版本回退到3.1来解决. 因此我们需要手动编辑虚拟机的XML配置文件(需要打开直接编辑的选项), 找到对应的条目改成如下配置:

```
...
<type arch="x86_64" machine="pc-q35-3.1">hvm</type>
...
```

![game_xml_config.jpg](https://raw.githubusercontent.com/SilenWang/Gallary/master/game_xml_config.jpg)

终于...可以坐下来好好玩游戏了~~~~~

后记: 一篇9.8就开始写的博文硬生生被我鸽到了快10.8......猛汉王果然魅力巨大....