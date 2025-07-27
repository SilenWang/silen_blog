---
title: 再次尝试编译openfyde（下)
categories: Script
date: 2025-06-19 23:09:23
tags: ['硬件兼容', '编译', '镜像生成', 'Rockchip', 'fydeos', 'openfyde']
---

换了主板和CPU后，一切都顺畅了... 也不知道是原来主板还是CPU的问题，编译几十分钟就会报错，然后再进行编译又不会在上次的地方报错。着还是我人生第一次真遇到硬件不兼容这种问题... 再次活久见。

于是，编译openfyde的下篇这就来了，可惜在被破硬件折磨的过程中，r132-dev版本的 prebuilt image已经更新了... 没赶上热点呀...

<!-- more -->

其实我上次已经进行到最后几步了，只是在编译chorme这个庞然大物的时候因为我硬件的问题被卡住了。在完成编译后，就可以生成直接刷写到u盘的镜像了。

只是这一次我还真解决了一个上次的问题，即得到了`chromiumos.bin`之后，如何获得能够用Rockchip镜像工具刷写的镜像。按照起步走教程得到的是bin镜像，而官方的预编译都是img，明显格式不同。同时用Rockchip工具对官方的openfyde prebuilt镜像解包的话，似乎会得到一些不包含在bin中的内容。

其实官方是公布了工具的，它就是[openFyde/rk3588-image-maker](https://github.com/openFyde/rk3588-image-maker)，只是上次作为一个纯小白，都看不明白编译过程中到底在干些啥，也就更看不明白openfyde下的项目都是些啥了。

生成所需镜像的步骤很简单，还是先克隆这个项目，然后进根目录，跟着这个项目的说明运行下面的命令就行：

```bash
# 从bin中挂载需要的内容
./map_chromiumos_image.sh /PATH/TO/CHROMIUMOS/IMAGE.bin --board fydetab_duo
# 生成 update.img
./rk3588-mkupdate.sh
```

须注意，在我这次编译的r132版本上，项目下的 `Image/parameter.txt`文件需要修改，即`0xa0006d@0x00b02040(STATE)`要改成`0xa0006d@0x00b02040(STATE)`

之后将生成的image拷贝下来，用rockchip官方工具刷写就好。这次编译的镜像触屏、蓝牙、Linux容器都正常，应该是没什么问题了。

接下来就可以试着往里头塞东西啦！
