---
title: Using minicap to Capture Screenshots from a Phone
categories: Others
date: 2018-08-28 12:05:17
tags: ['android', 'adb', 'minicap']
---

Practicing automated phone operations inevitably involves capturing the current image of the phone. The phone's built-in screenshot is too slow, so replacing it with minicap.

<!-- more -->

# Software Compilation

```bash
# Clone the project code and initialize
git clone https://github.com/openstf/minicap
cd minicap/
git submodule init
git submodule update
# Compile, ndk is required (archlinuxcn source has it, just install directly)
sudo pacman -S android-ndk
# Compile, after installation, it won't be added to bin automatically, so run with absolute path
/opt/android-ndk/ndk-build
```

# Software Deployment to Phone

```bash
# To use adb connect to the device, if not installed first install adb
# manjaro community has it, I remember the official version is自带的, now using the community edition may be due to minimalization considerations not default installing
sudo pacman -S android-tools

# Turn on USB debugging on the phone and connect to the computer, if device is not registered in phone select allow current device debugging
adb devices

# Get phone basic information to push corresponding version of minicap to the phone
ABI=$(adb shell getprop ro.product.cpu.abi | tr -d '\r')
SDK=$(adb shell getprop ro.build.version.sdk | tr -d '\r')
adb push libs/$ABI/minicap /data/local/tmp/
adb push jni/minicap-shared/aosp/libs/android-$SDK/$ABI/minicap.so /data/local/tmp/
adb shell chmod 777 /data/local/tmp/minicap
# Test if it can run
adb shell LD_LIBRARY_PATH=/data/local/tmp /data/local/tmp/minicap -P 1080x1920@1080x1920/0 –t
# If running successfully, the last will display 'OK'
```

# Screenshot

```bash
# The environment variable setting is necessary because it needs to specify the library path, -P后面的参数是截图参数, 详情可参见minicap的项目网页, -s but表截图并输出到
# 标准输出, 所以重定向就好
adb shell LD_LIBRARY_PATH=/data/local/tmp /data/local/tmp/minicap -P 1080x1920@1080x1920/0 –s > /sdcard/minicap/test.jpg
```

# To be completed
- Wrap the above commands into a script, and set up wifi connection information at the same time, so you can control multiple devices simultaneously.
```