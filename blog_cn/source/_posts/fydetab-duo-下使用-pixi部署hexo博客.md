---
title: fydetab duo 下使用 pixi部署hexo博客
categories: Script
date: 2024-08-04 13:04:05
tags: ['FydeOS', 'Fydetab Duo', 'arm', 'Linux环境配置', 'hexo', 'pixi']
---

<!-- 摘要部分 -->
<!-- more -->


Fydetab duo 也拿到了一个多月了，目前还在持续想办法把设备用起来的阶段。最起码，我希望能用这个设备来抽空谢谢博客。因此有了下面的折腾。

<!-- 摘要部分 -->
<!-- more -->

老实说，fydeos作为中国版本的chromeos，推荐的使用linux方法本应该是直接使用linux子系统的。但是在这一个月内，linux子系统似乎还存在一些影响使用的问题，比如唤醒后会莫名其妙的连不了网、子系统内鼠标显示不正常、子系统启动相对慢等比较影响使用的小毛病。加之fydeos团队对宿主系统的限制暂时还是比chromeos宽松的（至少使用sudo不需要切tty），因此我还是尝试在宿主系统部署我需要的工具。


## 安装pixi

- 安装`pixi`需要进行额外的环境变量配置，因为chromium os体系下，不是所有地方都可以存放并运行可执行文件，只有`/usr/local/share`下可以，因此需要在环境变量文件中增加下面的内容，然后再运行`curl -fsSL https://pixi.sh/install.sh | bash`

```txt
export PATH=/usr/local/share/pixi/bin:$PATH
export PIXI_HOME=/usr/local/share/pixi
```

## 安装code-server和git

原本首选的是微软官方给的`vscode-cli`工具，但可能是因为fydeos的文件结构和环境配置跟一般linux发行版还是有比较大的差异，因此实际上连接到服务时，会一直卡在下载vscode服务端的步骤，无法实际使用，因此暂时使用pixi安装`code-server`

另外fydeos的宿主系统已经自带了

## arm平台特殊配置

- fydeos tab是arm架构的设备，所以`pixi`的配置中需要增加`linux-aarch64`平台，同时实际测试下来yarn安装以来会出现无法识别子环境的python问题，因此fydetab duo不使用yarn，用回npm，同时由于fydeos中默认并没有`make`、`gcc`、`gxx`这些，因此也要将这些内容补充到`linux-aarch64`的依赖中
- fydeos的系统信息文件与常规linux有小差别，这会导致hexo的模块解析系统信息时错误，导致`hexo -v`和会调用它的`hexo g`无法运行，因此需要在部署任务中增加修复的语句

```toml
[target.linux-aarch64.tasks]
build_blog = {cmd = "npm install -g hexo-cli && npm install --python=$PWD/.pixi/envs/default/bin/python2 && sed -i '34s/distro\\[1\\]/distro/' $PWD/.pixi/envs/default/lib/node_modules/hexo-cli/dist/console/version.js", cwd = "."}
```

## 克隆代码并部署

要额外进行的配置就是这些了，把项目克隆下来直接部署，就可以继续写博客了

```bash
git clone https://github.com/SilenWang/silen_blog
pixi run build_blog
```
