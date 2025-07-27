---
title: 给Github首页做点面子工作
categories: Script
date: 2023-01-27 22:15:30
tags: ['Github', '个人主页', '美化', 'Github Stats', 'Shields.io', 'README']
---

给博客升级了, 那也要给Github首页也做点面子工程...
<!-- 摘要部分 -->
<!-- more -->

Github目前可以创建一个特殊的Repo, 该Repo和用户的名称一致, 比如我的就是`SilenWang/SilenWang`, 在该Repo下的`README.md`可以显示在用户的Gtihub首页上, 可以理解给了一个可以自定义的看板

目前我在上面放了一个Github的统计和我常用的工具图标

## Github统计

- 这个很简单, 是用一个API实现的:

```html
<img width="500px"  alt="GitHub Stats" src="https://github-readme-stats.vercel.app/api?username=SilenWang&count_private=true&show_icons=true"/>
```

## 工具图标

- 这个使用了[Shilds.io](https://shields.io/)这个工具, 可以通过get请求从网站API处请求回各种带图标的标签, 比如Jupyter的标签就用:

```html
<img alt="Jupyter" src="https://img.shields.io/badge/-Jupyter-f37524?style=flat-square&logo=Jupyter&logoColor=white" />
```

- 我把我用过的... 能查到有图标的都挖出来了...

- 可用的图标很丰富, 但也是不是什么都有, 目前已有的可以参见`simple-icons`下的[这个文件](https://github.com/simple-icons/simple-icons/blob/develop/slugs.md)

## 渲染效果

{% image https://raw.githubusercontent.com/silenwang/Gallary/master/2023/01/upgit_github_front_20230127_1674830298.png %}

## 后记

这边装饰好后, 我尝试把同样的标签用在`关于我`页面了, 然而Hexo渲染HTML Tag的时候似乎会强制加换行符, 导致所有标签拍成一列, 巨{% ruby 壮观|难看 %}... 我找了半小时没找到什么方法解决... 于是就摆烂用表格把每个标签框起来了... 以后找到如何解决再说吧...
