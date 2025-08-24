---
title: Compiling ARM version of Theia-IDE using GitHub Actions
categories: Others
date: 2025-06-16 23:20:58
tags: ['github', 'theia']

Since my first encounter with ChromeOS/FydeOS, I've been trying various VSCode-like editors. Recently, after learning about Huawei's CodeArt and its upstream project Theia, I started tinkering again. 
Unfortunately, my current device is an 8GB FydeOS, not the 16GB Manjaro or PixelBook 2017 I used before. The available memory in the Linux container is very limited - I can't even compile an ARM version of Theia-IDE browser edition... So... I had to "borrow" GitHub's resources again.

<!-- more -->

This combines what I mentioned before about {% post_link 白嫖codespace用来写hexo博客 [Using Codespace to write blogs] %} and {% post_link github-actions使用 [Using GitHub Actions] %}.

First, create a new project on GitHub, then click `Code` in the upper right corner to create a Codespace for this project, which allows you to use VSCode directly in the browser for coding and saving.
I directly referenced my previous [certimate_win7](https://github.com/SilenWang/certimate_win7) project to create a GitHub workflow, making sure to use GitHub's Arm Runner.

[The original Theia-IDE project](https://github.com/eclipse-theia/theia-ide) doesn't have official releases yet, so we can't download pre-built releases by tags. Instead, we need to use checkout to get specific code versions.

Then follow the instructions in the project's [browser.Dockerfile](https://github.com/eclipse-theia/theia-ide/blob/master/browser.Dockerfile) to set up the compilation steps. After that, package the entire code directory and upload it to my project's [Releases](https://github.com/SilenWang/theia-ide-browser-build/releases).

After downloading the compiled files, install system-level dependencies as per the official instructions, then run `yarn browser start` to launch it.

Since I'm using FydeOS, I compiled the browser version directly. If you need the desktop version, just follow the official instructions to generate the ARM version.
