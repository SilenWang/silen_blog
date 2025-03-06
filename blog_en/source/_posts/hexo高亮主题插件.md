---
title: hexo syntax highlighting plugin
categories: Others
date: 2018-09-05 21:53:50
tags: ['hexo', 'syntax']
---

I found that the default theme of Hexo is not very appealing, so I wanted to find a bright one and tried using a plugin.

<!-- more -->
The plugin can be found [here](https://github.com/ele828/hexo-prism-plugin). Installation is quite simple; just run `npm i -S hexo-prism-plugin`. Then, open the blog's configuration file `_config.yml`, disable the built-in highlighting, and add the settings for this plugin.

```yaml
highlight:
  enable: false
prism_plugin:
  mode: 'preprocess'    # realtime/preprocess
  theme: 'default'
  line_number: false    # default false
  custom_css: 'path/to/your/custom.css'     # optional
```

After that, update as usual. This plugin supports multiple highlighting modes and custom CSS. The supported languages are listed on the above link.
```