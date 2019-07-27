---
title: Git使用笔记
categories: Others
date: 2018-09-28 12:34:12
tags: ['git']
---

很早之前就看过Git的教程, 一直没实际使用过, 现在用起来, 感觉自己有那么点像程序员了呢~
<!-- more -->

Git 是一种分布式的版本管理程序, 方便多人协作/管理一个项目的代码. 不过对我来说, 因为绝大多数时候是单人使用, 所以对我来说更多的是在不同设备上同步我的代码了.

# 操作方法记录

- 用户初始化, 新设备上第一次使用时需要

```bash
git config --global user.name "Silen Wang" # 设置用户名称
git config --global user.email "mymail@gmail.com" # 设置用户邮箱
git config --global core.editor 'vim' # 设置提价commit时使用的
```

- 克隆一个项目的所有代码:

```bash
git clone PROJ_URL
```

- 拉取更新

```bash
git pull # 适用项目来自克隆的情况, 如果有多个来源好像要指定从哪个来源拉
git pull origin master # 指定拉取origin远程库的master分支上的修改
```

- 推送更新

```bash
git push
git push origin master # 指定推送到origin的master分支, 如果不指定好像会调用某种默认情况? 具体不是很清楚
```

- 远程管理

```bash
git remote # 查看远程信息
git remote add REMOTE_URL # 添加信息的远程库
```

- 修改提交, Git在对本地文件作出修改后, 修改是存在暂存区域的, 要提交这些更改, 首先得确认提交修改

```bash
git add FILE # 新增跟踪的文件/或将对文件的修改变成准备提交的状态
git commit # 提交所有未提交的修改, 这里会要求给一个简单的描述, 说明修改了什么东西, 通过这个描述可以回溯自己之前做了什么
git commit -a # 直接提交所有未提交修改, 并加上相同的描述
```

- 分支操作, 分支的概念允许创造当前项目的一个分支(就是一个备份), 然后在这个分支进行修改, 原分支(比如master)不受影响, 等修改的分支测试无问题后再将修改合并到需要修改的分支即可

```bash
git chcekout -b dev # 创建并切换到dev分支
git chcekout dev # 切换到dev分支
git merge dev # 将dev分支的改动合并到当前分支
```