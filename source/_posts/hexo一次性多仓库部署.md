---
title: 水一发:hexo一次性多仓库部署
categories: Others
date: 2020-08-02 01:19:38
tags: ['hexo', 'gitee']
---

最近github的访问越来越不行了, 有时候想查自己博客里的博文都费劲...所以在国内的git托管服务商gitee也托管了一份, 但是问题就来了...如何一次更新多地部署呢?
<!-- 摘要部分 -->
<!-- more -->

首先要设置在gitee上的仓库, 基本方式和github一致啦, 唯一需要注意的是gitee的默认展示项目和github不一样, github用`USERNAME.github.io`, 而gitee使用的是`USERNAME`.

由于gitee本身有fork github项目的功能, 因此在创建项目的时候选择fork原来在github上的博客仓库, 然后把仓库名改成`USERNAME`就可以(图里是迁移完成后举例的...所以显示已创建).

![gitee_proj](https://raw.githubusercontent.com/SilenWang/Gallary/master/gitee_proj.png)

之后再项目页面选服务菜单进行一下设置, 就能从`USERNAME.gitee.io`访问博客了.

![gitee_pages](https://raw.githubusercontent.com/SilenWang/Gallary/master/gitee_pages.jpg)

然后就是设置通过ssh等度gitee了, 这个也跟github的一致就不展开.

最后就是设置项目的`_config.yaml`文件, 在部署的部分如下配置, 这样就能部署的时候一次性传两次了.

```yaml
deploy:
  - type: git
    repo:
      git@github.com:SilenWang/silenwang.github.io.git
    branch: master
  - type: git
    repo:
      git@gitee.com:silenwang/silenwang.git
    branch: master
```

本来这样就该结束了...但是测试的时候发现另外一个问题...就是家里的网已经完全不能直连github, 导致部署不到github上了... 因此又搜索了下给git加代理的方式.

由于我配置里走的是ssh协议连接github, 经查需要通过ssh的配置文件添加一个`ProxyCommand`的设置, 如下:

```config
Host github.com
    User git
    HostName github.com
    ProxyCommand nc -x 127.0.0.1:1080 %h %p
```

添加之后再次尝试`hexo d`, 确实成功部署了.

以上, 水完~