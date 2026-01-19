---
title: 啊这...原来 devpod 已经被官方搁置了吗
categories: Others
date: 2026-01-17 02:54:00
tags: 
  - Devpod
  - 开源项目维护
  - 社区分支
  - 项目可持续性
  - Docker Compose
  - 开源软件
---

事情的起因是，使用 devpod 配置基于 docker-compose 的容器时，发现 devpod 似乎不能正确地调用 docker-compose 完成容器创建。结果不搜不知道，一搜... 啊？这项目居然被 loft-sh 放弃了吗？

<!-- more -->

根据 [issue 1946](https://github.com/loft-sh/devpod/issues/1946) 和 [issue 1915](https://github.com/loft-sh/devpod/issues/1915) 的描述，vCluster 是 loft-sh 商业化最成功、最有前景的项目，因此他们的所有精力都暂时放在上面，抽不出时间来维护社区项目了。

这个情况感觉跟 Fyde 团队类似，之前 Discord 上他们也说过，会优先倾向于能让公司继续运营的工作...

不过好在，一位[社区大神](https://github.com/skevetter)创建了 Devpod 的社区分支，并且以惊人的速度持续更新（目前已经 0.9.* 了）！

使用 skevetter 的 Fork，可以正常地使用 docker-compose，算是不幸中的万幸吧。

不过到目前为止，绝大部分的代码都是 skevetter 一人贡献的，这其实不是太健康。就如之前 onelist 那样，项目非常有名，但是几乎所有开发工作都集中在原作者一人身上，这样用爱发电的状态终究是无法持续的，一定程度上导致了最后卖项目跑路的问题。

希望 devpod 能像 gitea 那样，衍生自 gogs，但是最后发展出了完备的社区，让好项目能一直延续下去。
