---
title: 白嫖codespace用来写hexo博客
tags: ['云端开发', 'Hexo', '博客写作', 'Github', 'Codespace']
categories: Others
date: 2023-07-23 17:10:23
---


最近出于.....emmmmmm心有不甘？我卖掉了之前拍的lenovo chrome duet2，但是转头又买了pixelbook 2017.....既然设备到手了，还是想要继续折腾一阵，力求找到用这台pixelbook 2017代替我的工作电脑进行一般办公的方式。然后找着找着，主要矛盾没解决，倒是诞生了几个实用的副产品.... 使用codespace写hexo博客是其一。

<!-- more -->

## 在哪都能写博客啦
Github codespace本来是云时代，github推出的一项纯云开发服务，这个服务可以为用户分配一个容器（可以自定义），然后直接将github上的项目代码拉入到容器中，给用户分配一个可以访问容器内容的vscode web客户端，这样就能直接在github提供的服务器上直接进行代码的编写和调试，也就是完全在云端进行开发了。项目的实现我猜和[coder](https://github.com/googlecreativelab/coder)类似。

这么好的东西，当然不会是免费的，目前github提供15G的免费空间，以及120h的核心时间（[详细见官方](https://docs.github.com/zh/billing/managing-billing-for-github-codespaces/about-billing-for-github-codespaces)）。这个存储和时间看上去不多，但是根据目前的使用来看，核心时间应该是不作计算就不算的，因此如果是非常轻量的项目，其实还是很充裕的。

这不就正适合用来写博客嘛！

我的Hexo博客代码量还没用的主题多，需要的计算也就是生成一下静态网页，要不了几秒。云上操作还能解决以前的博客只能在办公电脑上写的窘境（只有办公电脑我配置了环境，懒得再来一次），打开即用。而且试了一下，`hexo s`生成的静态网页预览甚至都可以正常访问！这样就能跟那种博客平台一样，有浏览器就能写啦～（虽然我鸽显然不是因为没有设备能写

![ui](https://raw.githubusercontent.com/SilenWang/Gallary/master/2023/07/upgit_20230723_1690102891.png)

## 使用上的问题

目前唯二的问题就是网络以及博客内的图片管理了。

### 1. 图片问题

原来写博客用过[picgo](https://picgo.github.io/PicGo-Doc/zh/guide/)和[upgit](https://github.com/pluveto/upgit)两个工具进行图片上传和自动链接获取，目前实测下来，原本最方便的vscode-picgo插件并不能使用，即使安装`xclip`后，也识别不到剪贴板里的图像内容，大概是因为它是在云端容器中，而有图像的剪贴板，是我本地的吧... `upgit`倒是有办法使用，将软件下载到chromeos下的linux容器中，做好配置就可以了。 

### 2. 网络

这个问题不用多说... 国外的服务想访问到都有一样的问题，目前来说编辑器的部分其实还好，应为不是即时响应的，因此码字相当流畅...就是不知道会不会写着写着发现保存不了。

但是命令行的部分还是感觉得到明显卡顿的，好在命令行也不多，就是个`hexo g`和`hexo d`，忍一忍也就过去了，大不了就现在不发布，改天再说嘛。
