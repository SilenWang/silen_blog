---
title: AI还是改变生活--让我能做以前很难起步的事情
categories: Coding
date: 2025-11-22 20:25:43
tags: ['hexo', 'Volantis', '博客', '国际化']
---

我还是挺喜欢现在的 hexo 主题的，从 2021 年使用它的前身`material-x`，到2023年升级到 `Volantis 5.x`，主题的功能已经完全满足我个人的需求了，以至于没有继续升级作者后续的6和7了。不过上半年利用 AI 对博客进行批量英文翻译后，一直有一些界面元素还是中文的问题。今年就不大动工，直接在 AI 指导下自己动手修一修。

<!-- more -->

## 问题背景

自从使用 AI 工具将博客内容批量翻译成英文后，我的博客就具备了中英双语版本。但是修改的时候出现俩问题：

1. 英文站中部分组件的切换效果无效，如关于我页面的选项卡，点击后无法切换
2. 主题界面上的许多元素仍然是中文的，页面最下侧的访问统计也总是会出现一段明显是加入多语言支持时忘记处理的固定中文内容。

在查看了主题的代码后，感觉前者是个配置问题， 后者则更多是 hexo 这个框架本身的问题导致。比如页面中很多标题和分类的标题， 都是直接编码在配置文件中的，Deepseek 说，hexo 下更简单的方式是直接准备两个配置文件。我现在的方案本来就是中文站一个文件夹，英文站另一个，上述问题在实际尝试后发现，都是修改`Volantis`代码中的配置文件就能解决的，所以我最后选择直接写俩补丁文件，放在我的博客项目下。

## 解决方案

### 1. 不能切换的问题

从浏览器的报错看，当我把英文站的`root`设置为`/en/`后，页面上的`app.js`这个文件的路径，被错误设置成了`/en/en/js/app.js`，这个文件应该是负责实现切换效果的，文件找不到，切换效果也就无了。

根据 Deepseek 的解释，发现修改 `_config.yml` 下的cdn设置就能解决问题，具体的diff情况如下

```diff
--- a/_config.yml
+++ b/_config.yml
@@ -50,7 +50,7 @@ cdn:
   # 以下配置可以覆盖 cdn.prefix,配置项的值可以为空，但是要使用CDN必须依据路径填写配置项的键
   set:
     js:
-      #app: /js/app.js
+      app: /js/app.js
     css:
       #style: /css/style.css # (异步加载样式)
 # 静态资源版本控制

```

### 2. 英文化不全的问题

这些问题更简单，把需要显示板块的`title`替换成英文，然后找到页面底部不用的部分，去掉就行了。

```diff
--- a/_config.yml
+++ b/_config.yml
@@ -266,12 +266,12 @@ article:
       # 文章创建日期
       date:
         icon: fa-solid fa-calendar-alt
-        title: '发布于：'
+        title: 'Published on: '
         format: 'll' # 日期格式 http://momentjs.com/docs/
       # 文章更新日期
       updated:
         icon: fa-solid fa-edit
-        title: '更新于：'
+        title: 'Updated on: '
         format: 'll' # 日期格式 http://momentjs.com/docs/
       # 文章分类
       category:
@@ -279,15 +279,15 @@ article:
       # 文章浏览计数
       counter:
         icon: fa-solid fa-eye
-        unit: '次浏览'
+        unit: 'views'
       # waline 文章评论数量
       walinecount:
         icon: fa-solid fa-comment-dots
-        desc: '条评论' # 条评论
+        desc: 'comments' # 条评论
       # artalk 文章评论数量
       artalkcount:
         icon: fa-solid fa-comment-dots
-        desc: '条评论' # 条评论
+        desc: 'comments' # 条评论
       # 文章字数和阅读时长
       wordcount:
         icon_wordcount: fa-solid fa-keyboard
@@ -523,7 +523,7 @@ sidebar:
       sticky: true
       header:
         icon: fa-solid fa-list
-        title: 本文目录
+        title: TOC
       list_number: false
       min_depth: 2
       max_depth: 5
@@ -540,7 +540,7 @@ sidebar:
       display: [desktop] # [desktop, mobile]
       header:
         icon: fa-solid fa-folder-open
-        title: 文章分类
+        title: Categories
         url: /blog/categories/
     # ---------------------------------------
     # tagcloud widget
@@ -549,7 +549,7 @@ sidebar:
       display: [desktop, mobile] # [desktop, mobile]
       header:
         icon: fa-solid fa-tags
-        title: 热门标签
+        title: Tags
         url: /blog/tags/
       min_font: 14
       max_font: 24
@@ -572,7 +572,7 @@ sidebar:
       display: [desktop]
       header:
         icon: fa-solid fa-award
-        title: 站点信息
+        title: Website ino
       type:
         article:
           enable: true
@@ -608,7 +608,7 @@ sidebar:
       display: [desktop, mobile]
       header:
         icon: fa-solid fa-clock WISTERIA
-        title: 最近更新
+        title: Laste Update
 ############################### Sidebar ############################### > end
 
 
@@ -647,8 +647,6 @@ site_footer:
   source: https://github.com/volantis-x/volantis-docs/
   # analytics using leancloud
   analytics: >
-    <span id="lc-sv">本站总访问量为 <span id='number'><i class="fa-solid fa-loader fa-spin fa-fw" aria-hidden="true"></i></span> 次</span>
-    <span id="lc-uv">访客数为 <span id='number'><i class="fa-solid fa-loader fa-spin fa-fw" aria-hidden="true"></i></span> 人</span>
   # site copyright
   copyright: '[Copyright © since 2017 XXX](/)'
   # You can add your own property here. (Support markdown, for example: br: '<br>')
```

### 3. 补丁生成

将上面修改保存后，`git commit`，然后就可以使用`git diff COMMIT1 COMMIT2 > PATCH_FILE`生成补丁。

由于我是克隆`volantis`原项目代码然后修改生成补丁，补丁文件中文件的路径是不对的，需要修改补丁文件中下面的路径为实际的文件路径。

```diff
--- a/_config.yml
+++ b/_config.yml
```

另外上方的`index 91bd9709..9304175e 100644`需要去掉，否则会因为commit核验补上而失败。

## 总结

今年二季度以来，因为部门逐渐只剩我一个人，公司所有需要自己维护的网站和系统都跑我手上了... 如果没有AI救我狗命，我断然是无法一边做着分析，一边还要维护这么些的前后端的，现在也不可能有那么点常识后，自己修博客的小Bug... 但就跟去年用AI的感受一样，我的上限决定AI的上限，我们有办法判断拿到方案对错的时候，就只能像追自己尾巴的小狗狗一样，原地转圈圈了...
