---
title: 换新主题&Hexo框架搭建博客记录
categories: Others
date: 2019-03-12 00:20:02
tags: "Hexo"
---

昨天因为markdown文章里面用了无法解析的符号, 导致我以为hexo挂了, 遂折腾换了新的主题...顺便又回顾了一下搭建过程
<!-- 摘要部分 -->
<!-- more -->


## 一般搭建过程
hexo基于Node.js, 在arch系系统下`sudo pacman -S nodejs npm`即可

安装后使用npm在全局环境中安装hexo和hexo-cli(hexo-cli必要, hexo是否必要线下次再测)

```bash
# -g代表全局模块, 安装后可使用hexo命令行工具
sudo npm install -g hexo
sudo npm install -g hexo-cli
```

接下来使用hexo创建新的博客项目, 一句`hexo init`就好, 不过要注意, 这一步需要在空目录进行, 否则会提示目录非空, 无法运行.

hexo创建的博客有个默认的主题, 不过一般我们不会去用它就是了...在网上搜索一个你喜欢的主题, 然后使用git克隆到项目的theme目录下即可

之后对项目下的`_config.yml`进行修改就可以修改博客设置了. 修改完毕后`hexo g & hexo s`可生成并本地预览博客. `hexo d`则可以发布/跟新博客网站, 不过这需要在配置文件内指定好使用的方式(我是git), 并按照提示使用npm装好相应的发布用模块(直接`npm install`, 这个不需要装到全局)

## 博客项目备份

在很多情况下, 都可能会涉及到在多个地方编辑博客内容的问题. 这个就很简单了....直接git整个项目就可以了, 唯一需要注意的是, 如果使用git安装主题, theme下每个文件夹本来就是一个git项目, 我个人觉得在建立整个项目的git时在`.gitignore`中设置主题文件夹内所有文件被忽略, 然后主题的`_config.yml`另外备份一下比价好

![new_blog](https://raw.githubusercontent.com/SilenWang/Gallary/master/new_blog.png)