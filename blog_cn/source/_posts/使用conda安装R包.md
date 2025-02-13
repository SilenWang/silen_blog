---
title: 使用conda安装R包及jupyter
categories: Others
date: 2018-08-28 22:37:57
tags: ['R', 'conda', 'jupyter']
---

由于工作需要, 得安装一篇文献中提及的一个包, 本两想着两条命令完事, 谁知道一大堆依赖要解决...然后linux下的R包都是要编译的, 这可是要了老命了...还好有conda这个东西...

<!-- more -->

说起来原来看公众号的时候还不知道为啥里面推荐的用conda来装需要的工具...这次算是体会了linux的依赖有多烦人...而且作为普通用户又没有完全权限...路被赌死了不少...像conda这种直接在自己山头自立为王的方式...虽然浪费空间, 可是非常又用啊!

# conda安装

- conda有两个...emmm我也不知道能不能说是发行版? 一个是ancoda, 另一个是minicoda, 可以理解是完整安装和最小安装的区别, ancoda自带一套异常完整的python环境, 连jupyter都打包好了, 而miniconda除了最基本的pyton, 好像没太多别的东西了, 因为我主要是用来装r的, 所以选择了[miniconda](https://conda.io/miniconda.html)
- 下载到的linux安装文件是一个异常大的shell脚本, 直接运行即可开始安装
- 安装完之后即可调用conda命令安装想要的东西, 比如我就可以开始安装r了


```bash
conda install r
```

- 调用conda安装的东西全部都由conda统一放置在conda的安装目录下, 里面的文件夹分布好像跟根目录的不是特别一样
- 在最初完成conda安装时, 如果不作特别制定, conda会改写原本的`~/.bash`文件, 在里面把conda的bin目录添加到`PATH`变量中, 但是优先级在系统默认的目录(`/bin`)之后, 因为我是要用conda的覆盖系统默认配置, 所以手动把这一行加到了`~/.bash_profile`中并把conda的bin目录调到了系统默认的前面
- conda在国内也是有镜像源的, 编辑`~/.condarc`(没有可以直接创建), 加入下列内容即可使用中科大的镜像(清华的看人说有问题, 就没用了). 当然这只是个示例, `pkg`下还有其他的源, 需要的可以都添加上


```yaml
channels:
  - https://mirrors.ustc.edu.cn/anaconda/pkgs/free/
  - defaults
show_channel_urls: true
```

- 之后进入R, 在里面再安装相关包就好. 这时安装的R包不会放在`~/R`下, 也是全部在conda的目录下面, 具体可参考安装时的输出.
- 由于R包都是要编译的, 我又是用的minicoda, 所以难免会有依赖不满足. 不满足的包可以使用`conda install`安装, 另外也可以上anconda的网站先所有有哪些版本, 再对应执行命令安装即可


本次安装在conda的助力下解决了绝大多数的依赖问题, 唯一不能解决的我在本地编译后上传了...花了两天总算是把包装好了...但愿用不要出什么问题

# jupyter安装与配置

R放在集群上, 调试起来多少有些不便, 毕竟我习惯了Rstudio或者jupyter notebook的那种编写边测试的方式. 正好conda也能装jupyter, 就搜索了一下远程访问的配置方式.

- 首先肯定使用conda安装jupyter了

```sh
conda install jupyter
```

- 完成后对jupyter进行配置


```sh
# 生成配置文件
jupyter notebook --generate-config
# 设定密码, 根据提示输入两次即可
jupyter notebook password
```

- 使用编辑器打开生成的配置文件(普通用户的话是, 具体看生成文件后的提示`~/.jupyter/jupyter_notebook_config.py`)
- 找到下面几项, 取消注释后更改相应数值, 其中密码见生成时提示的那个json文件, 里面会有一段长长的字符串, 全部复制贴过来就好

```txt
c.NotebookApp.ip='*'
c.NotebookApp.password = u'your_key_str'
c.NotebookApp.open_browser = False # 代表启动notebook服务时不打开浏览器并访问
c.NotebookApp.port = 8888 # 这个可以不指定, 会自动分配一个端口
```

- 然后就可以在能连接到集群ip的电脑上打开浏览器以`0.0.0.0:8888`的方式访问了, 输入密码即可使用