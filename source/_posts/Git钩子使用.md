---
title: Git钩子使用
categories: Script
date: 2019-08-17 19:33:27
tags: ['git', 'githook']
---

虽然我也用Git进行部分项目的代码管理, 但是其实很多Git的特性我都还不熟悉, 也没有用, 今天尝试了一下Githook的使用.
<!-- 摘要部分 -->
<!-- more -->
Hook, 直译就是钩子, 指的是预先定义好的一系列脚本. 这些脚本在我们使用特定的Git命令完成一些代码管理操作时自动执行, 以达到减少重复操作, 提高效率的目的.

比如我的博客源代码就是用Git管理的. 我现在写一篇博客并发布的步骤是:
1. `hexo new post`创建一篇新的博客文章的makrdown草稿
2. 打开新建的markdown文件, 编写内容
3. 编写完成后运行`hexo g`和`hexo d`将博客内容生成并推送到Gitpages的源上去
4. 运行`git commit`提交博客源代码的修改, 然后`git push`到源代码源上去

其中, 撰写博客的过程可能会有一些不固定的操作, 但是最后两部的内容其实是每次写博客一定会重复且不太会变动的步骤. 这些步骤就可以利用钩子来执行.

首先要来判断这些东西应该在什么时候执行: 我的博客源代码在`commit`之前, 我会首先确认内容. 所以理论上, 我使用任何一个`git commit`会触发的钩子都可以. 实际选用的钩子为`post-commit`, 这个钩子会在`commit`后执行, 所以我就把如下内容写到项目的`.git/hooks/post-commit`文件中(名字不要改, 特定的名称才会被执行).

```bash
#!/bin/bash
cd /home/silen/git_proj/silen_blog/
hexo g
hexo d
```

然后给这个文件加上可执行权限, 这样一来, 每次我对这个项目完成一次commit后, 就会调用这个脚本自动进行博客的生成和推送了.

当然现在这里举的例子可能没啥吸引力, 因为这两步其实也并不多...省不了太多功夫. 但是我现在的站点除了主站, 我还在子站里弄了个个人的速查手册. 这个手册使用mkdocs进行编写, 是另外一个独立项目. 如果想要达到两个项目其中一个更新就同步更新到博客的话, 钩子应该就很实用了.

最后需要注意的是, git钩子是分为客户端和服务端的, 我现在使用的是客户端的一种钩子. 客户端钩子的一个缺陷是无法在项目的不同副本间同步(不同副本是不同的客户端), 这个问题我之后再找方法解决好了...

以上~