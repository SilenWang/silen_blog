---
title: hexo代码高亮插件
categories: Others
date: 2018-09-05 21:53:50
tags: ['hexo', '代码高亮']
---

感觉hexo默认的主题不是太好看, 想找个亮色的, 所以找了个插件来试用.

<!-- more -->
插件的地址在[这里](https://github.com/ele828/hexo-prism-plugin), 安装相当简单, 直接`npm i -S hexo-prism-plugin`, 就行. 然后打开博客的配置文件`_config.yml`, 关闭自带的高亮, 在下面添加这个插件的设置字段就ok了.

```yaml
highlight:
  enable: false
prism_plugin:
  mode: 'preprocess'    # realtime/preprocess
  theme: 'default'
  line_number: false    # default false
  custom_css: 'path/to/your/custom.css'     # optional
```

之后像往常一样更新就行, 该插件支持多种高亮方式及自定义css, 支持的语言在上面的地址中也有.
