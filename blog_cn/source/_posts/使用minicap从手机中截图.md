---
title: 使用minicap从手机中截图
categories: Others
date: 2018-08-28 12:05:17
tags: ['安卓', '截图', 'minicap', 'adb', 'android']
---

练习手机自动操作不可避免要涉及到截取手机当前的图像, 手机adb自带的截图实在太慢, 所以安装minicap取代之.

<!-- more -->

# 软件编译

```bash
# 克隆项目代码并进行初始化
git clone https://github.com/openstf/minicap
cd minicap/
git submodule init
git submodule update
# 编译需要ndk(archlinuxcn源内有, 直接安装即可)
sudo pacman -S android-ndk
# 进行编译, 书面包安装之后不会自行添加到bin, 需要给绝对路径运行
/opt/android-ndk/ndk-build
```

# 软件部署到手机

```bash
# 要使用adb连接到设备, 如果没有安装adb首先安装
# manjaro community里面有, 我记得官方版本是自带的, 现在用的社区发行版可能处于精简考虑没有默认安装
sudo pacman -S android-tools

# 手机打开usb调试然后连接到电脑, device现实未注册的话在手机上选择允许当前设备的调试
adb devices

# 获取手机基本信息以推送对应版本的minicap到手机
ABI=$(adb shell getprop ro.product.cpu.abi | tr -d '\r')
SDK=$(adb shell getprop ro.build.version.sdk | tr -d '\r')
adb push libs/$ABI/minicap /data/local/tmp/
adb push jni/minicap-shared/aosp/libs/android-$SDK/$ABI/minicap.so /data/local/tmp/
adb shell chmod 777 /data/local/tmp/minicap
# 测试是否可运行
adb shell LD_LIBRARY_PATH=/data/local/tmp /data/local/tmp/minicap -P 1080x1920@1080x1920/0 –t
# 如果运行成功最后会显示'OK'
```

# 截图

```bash
# 前面的环境变量设定是必须的, 因为要指定运行库位置, -P后面是截图参数, 详情可参见minicap的项目网页, -s但表截图并输出到
# 标准输出, 所以重定向就好
adb shell LD_LIBRARY_PATH=/data/local/tmp /data/local/tmp/minicap -P 1080x1920@1080x1920/0 –s > /sdcard/minicap/test.jpg
```

# 待完成
- 把上述命令封装成脚本, 然后顺便设定wifi连接的信息, 这样可以实现同时控制多个设备
