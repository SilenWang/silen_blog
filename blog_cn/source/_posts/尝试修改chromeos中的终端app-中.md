---
title: 尝试修改chromeos中的终端app(中)
categories: Coding
date: 2025-09-05 01:53:10
tags: ['chromeos', 'terminal', 'libapp', '终端']
---

书接上回，在成功编译libapp项目下的代码后，接下来是两部分：修改libapp的代码，以及将修改应用到镜像中
<!-- more -->

## terminal app原理
上次我已经成功在浏览器中看到了`terminal app`，要修改代码，先要大致理解这个应用的原理。

在进行了一点尝试后，我对终端应的理解如下：

1. 终端用由几个页面组成，最外侧的界面实际上是 `terminal.html` 这个页面，页面调用，然后调用一系列项目内的js脚本加载页面内容
2. 建立ssh链接时，应用实际上是，调用chromeos下才可调用的本地API（类似electron？），打开一个新的窗口，然后显示`terminal_ssh.html` 这个页面，该页面会调用nassh目录中的js脚本来进行实际的终端渲染
3. 启动ssh连接，实际上是将ssh的参数传递给`terminal_ssh.html`页面，然后页面内的脚本拼凑出ssh命令，交由下游的api去建立ssh链接
4. 设置也是独立的页面，对应`terminal_settings.html`

## 修改libapp代码
了解了应用的工作原理，接下来就是尝试修改代码了。由于我并没有专门学习过html/js/css这些网站的技术栈，虽然大致能理解代码的逻辑，但是我并不打算自己从头设计并编写代码，主要还是参考现有的代码，组合出我需要的功能。

就我自己的使用体验来看，sftp挂载到本地的这一现有功能，逻辑上和我想要的快速端口转发有一定相似性：

1. sftp挂载需要先设定ssh部分的端口、用户等信息
2. 在已有的登录条目基础上，进行sftp挂载，点击挂载后，与ssh登录一样新跳出一个页面，但是不能操作，只提示挂载信息，和关闭窗口后挂载自动结束

我想要的功能就参照这个逻辑来增改：
1. 端口转发需要先设定ssh部分的端口、用户等信息
2. 在已有的登录条目基础上，端口转发，点击转发后，跳出窗口设置要转发哪个端口
3. 与ssh登录一样新跳出一个页面，但是不能操作，只提示转发的端口，和关闭窗口后挂载自动结束

[实际修改的文件](https://github.com/SilenWang/libapps/tree/feat/add-port-forward-button)主要是terminal个nassh目录下的js文件，主要是加入转发的按钮，然后是参考sftp的功能实现，改一个能执行端口转发的函数出来

## 将修改应用到镜像
在chromeos的源代码中，libapp并不是包名，在检索项目后，发现这些内容的编译实际上在`crosh-extention.ebuild`中，找到这个ebuild的目录，使用git生成补丁，放到这个目录中

```bash
git diff COMMIT_HASH_1 COMMIT_HASH_2 > patchfile.patch
```

然后还要进一步的修改ebuild文件，加入下面的内容后，构建时会自动调用补丁，进行代码修改

```ebuild
PATCHES=(
    "${FILESDIR}"/0001-forward.patch
)
```

## 重新编译包和镜像

这里我碰到了一个坑，如果只修改文件，直接按照之前的方法进行编译，其实是无法应用改动的。

根据[找到的这个内容](https://www.jianshu.com/p/6d8523b1f771)，需要使用`cros_sdk`中的`cros_workon`命令，标记对特定包进行修改，才会使用我修改后的ebuild文件来进行编译，而不是获取特定版本的源代码进行编译。

编译的过程{% post_link 编译openfyde [参考以前的步骤] %}就可以了。在执行`build-packages`时，工具会自动识别有更改的包，自动重新构建，然后将更新后的包放到chroot环境中去。

## 刷机以及启动

{% post_link 再次尝试编译openfyde-下 [刷机也没有什么变化] %}，就是上次随便改内核，导致测试用的工程机duo砖了，因此参考[fydetab duo wiki的说明](https://wiki.fydetabduo.com/unbrick_the_fydetab_duo)进maskrom模式来刷。


## 进行调试？

这里其实留下了一个问题，不论做什么开发，在进行实机测试和调试之前，理论上都应该先在测试或者虚拟环境中完成测试才是，否则要来回编译、刷机，会浪费大量时间。我想这也是各种SDK和Studio套件存在的意义。

但是，我目前还真不知道... 我这样修改应该要如何在刷机前调试？根据[找到的博客](https://www.owalle.com/2020/06/03/crosvm-chromevm/)，似乎是可以利用cros_sdk的虚拟机功能的，调用QEMU的kvm虚拟机，似乎可以对编译后的镜像进行调试。

下次就是把功能调试好，展示成果了！这次先展示一下刷机后界面的改变！

![terminal_in_openfyde.png](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/09/upgit_20250907_1757177833.png)
