---
title: Sprucing Up My GitHub Profile
categories: Script
date: 2023-01-27 22:15:30
tags: ['GitHub']
---

After upgrading my blog, I decided to do some cosmetic work on my GitHub profile...

<!-- Summary section -->
<!-- more -->

Currently, you can create a special repository that matches your username, like `SilenWang/SilenWang`. The `README.md` in this repository will be displayed on the user's GitHub homepage, essentially giving them a customizable dashboard.

I currently have my GitHub statistics and commonly used tool icons displayed there.

## GitHub Statistics

This is quite simple; it uses an API:

```html
<img width="500px" alt="GitHub Stats" src="https://github-readme-stats.vercel.app/api?username=SilenWang&count_private=true&show_icons=true"/>
```

## Tool Icons

I used [Shields.io](https://shields.io/) for this. It allows you to request various labeled icons from the website's API via a GET request. For example, the Jupyter label is:

```html
<img alt="Jupyter" src="https://img.shields.io/badge/-Jupyter-f37524?style=flat-square&logo=Jupyter&logoColor=white" />
```

I dug out all the icons I could find for tools I use...

The available icons are quite rich, but not everything is available. You can see the existing ones in [this file](https://github.com/simple-icons/simple-icons/blob/develop/slugs.md).

## Rendering Effect

{% image https://raw.githubusercontent.com/silenwang/Gallary/master/2023/01/upgit_github_front_20230127_1674830298.png %}

## Postscript

After decorating it, I tried to use the same labels on my "About Me" page. However, Hexo seems to force line breaks when rendering HTML tags, causing all labels to be in a single column, which is {% ruby ugly|unattractive %}... I spent half an hour trying to find a solution but couldn't figure it out. So, I just used tables to box each label for now. I'll look into how to fix this later.
