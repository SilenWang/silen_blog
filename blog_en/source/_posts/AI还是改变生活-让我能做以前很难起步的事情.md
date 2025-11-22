---
title: AI is Still Changing Lives - Enabling Me to Do Things That Were Once Hard to Start
categories: Coding
date: 2025-11-22 20:25:43
tags: ['hexo', 'Volantis', 'blog', 'internationalization']
---

I'm still quite fond of my current Hexo theme. From using its predecessor `material-x` in 2021 to upgrading to `Volantis 5.x` in 2023, the theme's functionality has completely met my personal needs, to the point where I haven't upgraded to the author's subsequent versions 6 and 7. However, after using AI to batch-translate my blog content into English earlier this year, there have been persistent issues with some interface elements remaining in Chinese. This year, instead of undertaking major renovations, I decided to fix these minor bugs myself under AI guidance.

<!-- more -->

## Problem Background

Since using AI tools to batch-translate my blog content into English, my blog has had both Chinese and English versions. But during the modification process, two issues emerged:

1. Some component switching effects were invalid in the English site, such as tabs on the "About Me" page that couldn't switch when clicked
2. Many elements on the theme interface remained in Chinese, and the visitor statistics at the bottom of the page always displayed fixed Chinese content that was clearly forgotten during the multilingual support implementation.

After examining the theme's code, I felt the former was a configuration issue, while the latter was more due to the Hexo framework itself. For example, many titles and category headings on the page were directly encoded in the configuration file. Deepseek suggested that a simpler approach in Hexo would be to prepare two configuration files directly. My current setup already uses separate folders for the Chinese and English sites, and after practical testing, I found that both issues could be resolved by modifying the configuration files in the `Volantis` code. Therefore, I ultimately chose to write two patch files and place them in my blog project.

## Solutions

### 1. The Switching Issue

From browser error messages, when I set the English site's `root` to `/en/`, the path for the `app.js` file on the page was incorrectly set to `/en/en/js/app.js`. This file should be responsible for implementing the switching effects. When the file couldn't be found, the switching effects were lost.

According to Deepseek's explanation, modifying the CDN settings in `_config.yml` could solve the problem. The specific diff is as follows:

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

### 2. Incomplete English Localization

These issues were simpler to fix - just replace the `title` of the display sections with English and remove the unused parts at the bottom of the page.

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

### 3. Patch Generation

After saving the above modifications, `git commit`, then you can use `git diff COMMIT1 COMMIT2 > PATCH_FILE` to generate a patch.

Since I cloned the original `volantis` project code and modified it to generate the patch, the file paths in the patch file are incorrect. You need to modify the following path in the patch file to the actual file path:

```diff
--- a/_config.yml
+++ b/_config.yml
```

Additionally, the `index 91bd9709..9304175e 100644` above needs to be removed, otherwise it will fail due to commit verification.

## Conclusion

Since the second quarter of this year, as my department gradually dwindled to just me, all the websites and systems that the company needs to maintain have fallen into my hands... If AI hadn't saved me, I absolutely wouldn't have been able to handle analysis work while maintaining all these front-end and back-end systems simultaneously. Nor would I have been able to fix these minor blog bugs with just a bit of common knowledge now... But just like my experience with AI last year, my upper limit determines AI's upper limit. When we can't judge whether the solutions we get are right or wrong, we can only spin in circles like a little dog chasing its own tail...
