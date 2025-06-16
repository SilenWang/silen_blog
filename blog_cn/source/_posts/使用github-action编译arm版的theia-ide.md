---
title: 使用github_action编译arm版的theia_ide
categories: Script
date: 2025-06-16 23:20:58
tags: ['github', 'theia']
---

自从第一次接触chromeos/fydeos开始，我就一直在尝试各种不同的 vscode(-like) 编辑器，最近因为看到华为的CodeArt，知道了它的上游项目theia，又开始了折腾。
奈何我手上的设备已经是只有8G内存的 fydeos 不是之前的16G Manjaro，也不是PixelBook 2017。在Linux容器下的可用内存十分受限，连编译个arm版的 theiaide-ide browser版都不得行... 那... 只有又来白嫖Github了

<!-- more -->

这次结合之前提到的 {% post_link 白嫖codespace用来写hexo博客 [白嫖codespace写博客] %} 和 {% post_link github-actions使用 [github action的使用] %} 。

首先在github上新建一个项目，然后找到右上角的`Code`创建该项目的codespace，即可在浏览器中使用vscode，直接进行代码编写和保存。
我直接参照之前已经写过的[certimate_win7](https://github.com/SilenWang/certimate_win7)项目，创建一个github workflow，注意Runner用Github提供的Arm Runner。

[Theia-IDE的原项目](https://github.com/eclipse-theia/theia-ide)还没有正式的Release，所以不能根据Tag下载打包好的Release，改成使用checkout来获得特定的代码。

然后按照项目下的[browser.Dockerfile文件内容](https://github.com/eclipse-theia/theia-ide/blob/master/browser.Dockerfile)，设置编译步骤，接着将整个代码目录打包，上传到我自己项目的[Release中](https://github.com/SilenWang/theia-ide-browser-build/releases)。

下载编译文件后，容器内按照官方的说明安装好系统级别的依赖，`yarn browser start`就可以启动了。

因为我是fydeos用，直接编译的browser版，如果需要桌面版本，按照官方的说明生成Arm版就好了。