---
title: Trying to Modify the Terminal App in ChromeOS (Part 1)
categories: Coding
date: 2025-08-21 00:13:10
tags: ['chromeos', 'terminal', 'libapp']
---

I've been using FydeOS/ChromeOS for about two years now. While the system provides a usable terminal app, it's honestly not that great. For example, during development, I often need to forward multiple ports. Although I can achieve this by entering SSH commands for port forwarding, this requires manually inputting quite a few parameters. Additionally, during port forwarding, I need to keep the SSH login window open. For someone like me who's particularly obsessive about minimizing the number of open windows, keeping three or four windows open that won't be used in the foreground is really uncomfortable... So I thought, can I do it myself, with the help of AI, modify the system's default terminal client, and add quick forwarding functionality like VSCode has?

<!-- more -->

## Cloning the Code
The first step is to find out if the project code can be obtained. DeepSeek, Kimi, and ChatGPT can all give me the correct answer: `https://chromium.googlesource.com/apps/libapps`. Just clone this project locally.

## Project Compilation
Almost every directory in the project comes with documentation, some detailed and some brief. After reading through it and combining it with AI's answers, the key folders are roughly as follows:

- hterm: A terminal emulator written in JS. All the content you see in the terminal after logging in is rendered by it
- nassh: The terminal emulator extension in the Chrome Web Store, combining content from `hterm` and `ssh_client`
- ssh_client: The part that communicates with OpenSSH
- terminal: The frontend of the terminal app we see

## Trying to Run Outside ChromeOS
To develop/modify an application, you first need to set up the corresponding development and debugging environment. This was the most challenging part for me. Although the project documentation is quite rich, it's not rich enough for a beginner like me to easily learn how to build a test version of the app that can be debugged. So with the help of AI, I found at least a solution that can show the interface:

- Use pixi to create a virtual environment (On my machine, this is necessary, otherwise the next step will report file copy permission errors)
```toml
[workspace]
authors = ["Sylens Wong <qiumin14@163.com>"]
channels = ["conda-forge"]
name = "fydeos_dev"
platforms = ["linux-64"]
version = "0.1.0"

[activation.env]
PATH = "$PIXI_PROJECT_ROOT/libapps/libdot/bin:$PATH"
```
- Enter the virtual environment with `pixi shell`, then run the project's built-in `kokoro/build` script to complete the build and testing of all sub-projects (Testing must be run completely, as the test part is directly written in)
- Enter the `terminal` directory, then copy some files needed by this frontend that aren't in this project
```bash
cp ../node_modules/xterm/css/xterm.css css/ # Missing CSS
cp ../nassh/js/* js/ # Missing JS scripts
cp -r ../nassh/_locales ../ # Missing international language files
```
- Run `npm run start` in the project root directory to start the HTTP service, then enter the following addresses in the browser to see the terminal's frontend page:
    + `http://localhost:8080/terminal/html/terminal.html`: Main page
    + `http://localhost:8080/terminal/html/terminal_settings.html`: Settings page

![terminal](https://raw.githubusercontent.com/SilenWang/Gallary/master/2025/08/upgit_20250821_1755708797.png)

Although some content on the page is usable, such as adding SSH configurations, the core functionality like actually making SSH connections doesn't work because this is just the frontend of the app. The actual connection/terminal rendering and other content aren't here.

However, this is at least a small step forward~
