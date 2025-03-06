---
title: chromebook入坑记?
categories: Others
date: 2022-10-16 21:44:36
tags: ['chromeos', 'chromebook', 'fydeos']
---

I have always been interested in the status of domestic Linux systems (although I haven't used any long-term except for an ancient version of Linux Deepin), so when Flint OS was still called FydeOS, I paid attention to this Chromium OS-based special edition. This also led me to understand why Chrome OS, which has a mysterious high share overseas in the past, is popular.
<!-- 摘要部分 -->
<!-- more -->

Chrome OS is developed by Google and its features include being lightweight, simple, and heavily relying on cloud services. In other words, it's an operating system designed for low-performance netbooks, with the idea that most users don't need strong performance to run a lot of professional local software; they just need a device that can access the required services through the web. Therefore, Chrome OS was originally designed not to run any services locally, and all functions needed to be accessed through website services. The browser Pro Max Plus system (then received a lot of complaints).

From its development, it seems this idea might not be entirely wrong, but... with the addition of Linux subsystems and Android subsystems, completely cloud-based systems are still too early.

However, this doesn't mean it's too early for my work...

My work involves handling large amounts of gene sequencing data, which is too huge to process using mobile devices yet. Therefore, there has always been a need in my work to separate actual computational tasks from other work. On the other hand, many important tools I use, such as Jupyter, RStudio, and VSCode, already have mature server versions. After setting up corresponding services, I started thinking about whether it would be possible to complete daily work requirements using a Chromebook or a device similar to Chrome OS (as mentioned earlier).

So, in the following days, I first tried FydeOS(x86), then bought Lenovo Duet, and finally entered the Chromebook pit... Previously, I had tried setting up self-hosted cloud services to prepare for working on Chromebooks only. However, like the countless pitfalls I've encountered over 10 years using Linux distributions and even experiencing life-threatening losses before being able to use them smoothly, there are still many issues with using such a minimalist system for work at this stage... I'll list these problems one by one later...
