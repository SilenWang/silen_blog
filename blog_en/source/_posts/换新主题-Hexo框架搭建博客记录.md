---
title: Switching Themes & Hexo Framework Blog Setup Notes
categories: Others
date: 2019-03-12 00:20:02
tags: ["hexo"]
---

Yesterday, due to some unparsable symbols in my markdown post, I thought Hexo had crashed, so I took the opportunity to switch to a new theme...and reviewed the setup process again.

<!-- more -->


## General Setup Process
Hexo is based on Node.js. On Arch-based systems, simply run `sudo pacman -S nodejs npm`

After installation, use npm to install hexo and hexo-cli globally (hexo-cli is essential, whether hexo itself is necessary will be tested next time)

```bash
# -g stands for global module, after installation you can use hexo command line tools
sudo npm install -g hexo
sudo npm install -g hexo-cli
```

Next, use hexo to create a new blog project with just `hexo init`. Note that this needs to be done in an empty directory, otherwise it will prompt that the directory is not empty and cannot proceed.

The blog created by hexo comes with a default theme, but generally we don't use it...Search online for a theme you like, then clone it to the project's theme directory using git.

After that, modify the `_config.yml` file in the project to change blog settings. After modifications, `hexo g & hexo s` can generate and preview the blog locally. `hexo d` can publish/update the blog website, but this requires specifying the publishing method in the configuration file (I use git), and installing the corresponding publishing module as prompted (just `npm install`, this doesn't need to be installed globally)

## Blog Project Backup

In many cases, you may need to edit blog content in multiple locations. This is quite simple...just git the entire project. The only thing to note is that if you use git to install themes, each folder under theme is originally a git project. Personally, I think it's better to set all files in the theme folder to be ignored in `.gitignore` when establishing git for the entire project, and then separately backup the theme's `_config.yml`.

![new_blog](https://raw.githubusercontent.com/SilenWang/Gallary/master/new_blog.png)