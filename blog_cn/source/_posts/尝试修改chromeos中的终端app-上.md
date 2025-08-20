---
title: 尝试修改chromeos中的终端app(上)
categories: Coding
date: 2025-08-21 00:13:10
tags: ['chromeos', 'terminal', 'libapp', '终端']
---

我是用fydeos / chromeos也有两年了, 系统中虽然提供了能用的终端app, 但是说真的, 还不是那么好用, 比如我在开发时, 常常会需要进行多个端口的转发. 虽然我可以通过输入端口转发的ssh命令来达成目的,但是这样一来会需要手动输入比较多的参数, 二来端口转发期间会一直需要保持ssh登录后的窗口, 对于我这种特别喜欢降低开启窗口数目的强迫症来说, 开着三四个不会前台使用的窗口真的很难受... 所以我就想, 我能不能自己动手, 在AI帮助下, 修改系统自带的默认终端客户端, 给他加上vscode那样的快捷转发功能呢?

<!-- more -->

## 克隆代码
百尺竿头第一步, 先得知道项目代码是否能获取, 这个deepseek, kimi 和 chatGPT 都能给我正确答案: `https://chromium.googlesource.com/apps/libapps`, 直接克隆这个项目到本地就行.


## 项目编译
项目内几乎每个目录都自带说明, 有详细有简略的, 阅读一遍, 结合AI的回答可知, 关键的文件夹大致如下:

- hterm: js编写的终端模拟器, 登录后看到的终端所有内容都是它渲染的
- nassh: chorme商店的终端模拟器插件, 结合了`hterm`和`ssh_client`的内容
- ssh_client: 与openssh通信的部分
- terminal: 我们看到的terminal app前端


## 在chromeos外尝试运行
要开发/修改一个应用, 首先需要部署相应的开发及调试环境, 这个是这次最为难我的地方了, 虽然项目的文档内容比较丰富, 但是仍然不至于丰富到让小白如我能轻易学会怎么构建一个可以进行调试的测试版本app, 所以借助AI, 找到了一个至少能看到界面的方案:

- 使用pixi创建一个虚拟环境(在我的机器上必须这么做, 否则下一步会提示文件无权限复制的错误)
```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda-forge"]
name = "fydeos_dev"
platforms = ["linux-64"]
version = "0.1.0"

[activation.env]
PATH = "$PIXI_PROJECT_ROOT/libapps/libdot/bin:$PATH"
```
- `pixi shell`进入虚拟环境, 然后运行项目自带的`kokoro/build`脚本完成所有子项目的构建和测试(必须要跑完测试, 测试部分直接写进去了)
- 进入`terminal`目录, 然后将一些这个前端需要, 但是并不在这个项目内的文件复制过来
```bash
cp ../node_modules/xterm/css/xterm.css css/ # 缺少的css
cp ../nassh/js/* js/ # 缺少的js脚本
cp -r ../nassh/_locales ../ # 缺少的国际语言文件
```
- 在项目根目录运行`npm run start`启动http服务, 然后浏览器输入如下地址, 就能看到terminal的前端页面了:
    + `http://localhost:8080/terminal/html/terminal.html`: 主页面
    + `http://localhost:8080/terminal/html/terminal_settings.html`: 配置页面

![terminal](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/08/upgit_20250821_1755708797.png)

虽然页面上的一些内容是可以用的, 比如添加ssh配置, 但是基本的核心功能, 如实际进行ssh连接, 都是无法使用的, 因为这只是这个app的前端, 实际的连接/终端渲染等内容都不在这里.

不过, 总算是往前进了一小步~