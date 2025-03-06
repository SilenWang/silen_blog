---
title: bash-it
categories: Script
date: 2018-10-29 21:17:23
tags: ['bash', 'plugin']
---

bash-it is a collection of bash plugins inspired by oh-my-zsh. Installing bash-it not only makes the bash interface more stylish but also makes the command line more convenient.

<!-- more -->

## Installation

According to the installation tutorial on the bash-it homepage, it's very easy to install.

```bash
git clone --depth=1 https://github.com/Bash-it/bash-it.git ~/.bash_it
~/.bash_it/install.sh
```

However, when choosing, you need to pay attention that it will ask whether to retain the original `.bashrc` and `.bash_profile`. Of course, if you don't retain these two files, they won't be deleted; only backed up. If there are any problems after installation, you can manually switch back.

## Usage

`bash-it -h` can display all available commands, and `bash-it show plugins` can show all available plugins. Currently, I have only used `base` and `git`, but I will gradually learn about the others later.

![display](https://raw.githubusercontent.com/SilenWang/Gallary/master/bash_it.png)
