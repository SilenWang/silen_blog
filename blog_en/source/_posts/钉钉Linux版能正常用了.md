---
title: DingTalk for Linux Just Works Now
categories: Daily
date: 2022-02-22 22:22:22
tags: ['DingTalk']
---

Although DingTalk is somewhat of a tool that exploits people, it's still necessary for work. It's already better to have one less annoyance choice (speaking of which, Tencent). After all, being exploited is bad enough, but being both exploited and annoyed... being able to eliminate one annoyance is already an improvement.

<!-- more -->

Thanks to [FlyInWind](https://aur.archlinux.org/packages/dingtalk-bin) on AUR providing the solution, the Linux version of DingTalk can now be used normally on Arch-based systems. Currently, all components I need are working fine (chat, approval in the workbench, documents, project management, online document analysis). Video conferencing should be tried this week; if it also works... then everything will be perfect.

Currently, the developers have likely incorporated FlyInWind's solution directly into the build files. The DingTalk installed via AUR automatically deletes `libgtk-x11-2.0.so` and other two files. As long as you add the following code to the startup script according to dbh625's instructions, it will work perfectly.

```bash
export XMODIFIERS="@im=fcitx"
export QT_IM_MODULE="fcitx"
export QT_QPA_PLATFORM=xcb
```

Looking forward to the day when all commonly used software has a fully functional Linux version.
