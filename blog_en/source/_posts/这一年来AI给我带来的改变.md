---
title: How AI Has Changed My Life This Year
categories: Others
date: 2026-01-06 22:30:10
tags:
  - AI
  - Automation
  - DeepSeek
---

From late 2022 to early 2023, ChatGPT exploded in popularity, and I started using LLMs to assist with script writing. Although it was quite useful, the price was relatively high, and payment was always an issue, so its application was limited to simple coding problems. 

From late 2024 to early 2025, DeepSeek became a hit. While its answers often weren't entirely satisfactory, and the addition of a "thinking mode" made the response time a bit slower, it was really cheap! A year later, DeepSeek remains one of the most affordable models in terms of tokens while still delivering decent output. I paid 50 RMB at the beginning of the last year, and now, almost exactly a year later, I still have 27 RMB left... As a result, I've become much bolder in applying LLMs to various other areas.

<!-- more -->

## From "Luxury" to "Daily Necessity"

If using AI in the previous two years felt like an occasional luxury, this year it has completely transformed into a daily necessity in my workflow. DeepSeek's low cost allows me to use it without hesitation for all sorts of tasks: from simple command queries, to entire script refactoring, to exploring completely unfamiliar domains. While it may not have achieved a leap in efficiency, it has genuinely expanded my capabilities, enabling me to tackle work that I previously only had conceptual ideas about but never practiced (and then became an even more diligent workhorse).

## Successful Use Cases This Year

### Batch Translation of Blog Posts

Using the command-line AI tool [Aider](https://aider.chat/), I wrote a simple shell script that automatically reads each article from the Chinese blog directory, translates it into English, and saves it to the corresponding English directory.

The script completed the translation in just an hour and a half, but I spent a whole week of evenings manually proofreading... This is also why, despite having the idea to translate my blog into English in the past, I never followed through—without a script tool, copying and pasting article by article would have been exhausting. Additionally, that week of proofreading English translations reminded me of the time I spent preparing English clinical trial materials... Even though I've now read thousands of English professional papers and materials... I still find the process quite torturous... I guess that's just a lack of language talent...

### Contributing Chinese Translations to the FydeTab Duo Wiki

Besides my own blog, I also contributed Chinese translations to the FydeTab Duo Wiki several times. Again, I used Aider for the initial translation, then performed a round of proofreading myself. This also helped me learn about the FydeTab Duo Wiki, as I hadn't previously used the Wiki framework it use.

### Drafting Blog Posts

Starting in the second half of this year, the number of blog posts I uploaded began to increase exponentially, thanks to the combination of Aider + DeepSeek. Since I write blogs mainly for my own records, without much stylistic requirement—the purpose of documentation far outweighs that of presentation—I can just outline the points, throw the relevant code or Notebook at them, generate a draft, and then edit it. This really saves time. Previously, I could spend an entire evening on a single blog post... Now, as long as I have a topic, writing five posts in one evening isn't a problem... Combined with the earlier translation work, simultaneous Chinese and English publication is also feasible.

### Batch Generation of API Tests

Faced with Swagger description files from multiple service groups, manually writing pytest test cases used to be a not‑difficult but quite time‑consuming task. Again, with the help of Aider, I only needed to provide the Swagger description files, and it could generate a runnable basic test code framework within minutes. I just needed to fine‑tune on that basis to quickly improve test coverage.

When building a CI system for a system I took over in the first half of this year, the Aider + DeepSeek combination allowed me to complete the work much earlier than expected.

### Fixing Internationalization Issues in a Hexo Theme

While setting up the English site for my blog, I discovered that a small portion of the text in the Volantis theme was actually hard‑coded in configuration files or project code, and some interactive features failed under sub‑directory routing. By using DeepSeek to analyze browser console errors and the theme source code, I quickly pinpointed the root causes:
1.  Incorrect CDN path configuration caused `app.js` to fail loading; after fixing this, interactive functionality returned to normal.
2.  I replaced the fixed Chinese titles, descriptions, and other fields in the theme configuration file one by one with English equivalents.

In the end, I resolved the bilingual support issue by generating and applying patches, without directly modifying the theme source repository. This would have been very hard in the past, as I haven't really done full‑scale JS project development; reading code to find problems would have been quite laborious, but AI helped narrow down the inspection scope, making it possible to solve the problem in a short time.

### Configuring GitHub Actions for CI and CD

This is another good example. This year I gained experience with automation workflows on two platforms: Azure Pipeline and GitHub Actions. For the former, I wrote the configuration manually by following the documentation, which took about a full day, because the documentation was entirely in English, and the translation was rather poor (many of Microsoft's Chinese documents feel as if there are no Chinese users at all, so they are all machine‑translated without optimization).

GitHub Actions, on the other hand, was largely completed under AI guidance: generating a rough example and then me modifying it, basically taking 1–2 hours to finalize.

### Contributing ARM64 Packages to conda‑forge

When I discovered that `tree_sitter_languages` lacked a conda package for the ARM64 architecture, preventing `aider‑chat` from being installed on the Fydetab Duo (an ARM device), I decided to supplement it myself. Guided by AI, I:
1.  Learned about the feedstock mechanism of conda‑forge.
2.  Added `linux_aarch64` build support in the corresponding `conda‑forge.yml`.
3.  Followed the community process to submit a PR, which successfully passed CI tests.

## Failures or Partial Failures This Year

AI is obviously not omnipotent; its proposed solutions often require professional judgment, otherwise you can end up spinning your wheels, and many times I didn't even know whether I was spinning or just hadn't reached the end yet.

### Attempting to Modify Theia's Remote Module

I wanted to add remote development support to the browser‑based Theia IDE because using VSCode in the Linux subsystem on my FydeTab Duo was too laggy... Both Aider and OpenHands could provide code modification suggestions, but since I knew nothing about large Electron/Node.js projects, I couldn't judge whether those modifications were correct, nor how to debug them. In the end, this attempt came to nothing.

### Attempting to Update the FydeOS Kernel for Vulkan Support

To improve graphics performance, I considered upgrading the FydeTab Duo's kernel to a newer version in hopes of supporting Vulkan drivers. AI helped me find new kernel source code suitable for RK3588, but I couldn't understand the overlay replacement mechanism during compilation. Eventually, I took a crude replacement approach, which resulted in bricking the device.

### Attempting to Write a Custom Application Based on DevPod's Source Code

DevPod's functionality fits my requirement of creating a development environment independent of FydeOS's Linux subsystem, except that it must be based on Docker or Podman, and I currently cannot run Docker/Podman outside the subsystem. So I once tried to have AI directly reference and reuse DevPod's source code to write an application that connects via SSH, installs VSCode on a remote machine, and performs port forwarding. However, AI didn't reference or reuse any of DevPod's existing content; instead, it rewrote everything based on my description, and the resulting code didn't run.

### Attempting to Compile and Run Podman on FydeOS

I successfully used Pixi to manage dependencies and compiled an ARM64 version of Podman and its dependent libraries. However, at runtime, due to the SafeSetID security restriction in the FydeOS/ChromeOS kernel, user namespace mapping couldn't be completed, causing rootless Podman to ultimately fail to run. AI explained the origin of this issue to me, but to this day I'm not sure whether it was correct or not, nor how I should solve the problem.

## Summary

Looking back over this year, AI has become thoroughly integrated into my daily work and study. It has helped me complete a huge amount of repetitive labor, solved countless specific problems, and even given me the courage to venture into areas I previously wouldn't have dared to touch casually. But at the same time, I've become more acutely aware that the effectiveness of AI depends heavily on the user's knowledge base and judgment.

- **Within an existing knowledge framework**, AI can bring a certain efficiency boost. While it's absolutely not to the level of being able to be a hands‑off manager as some vendors claim, there is definitely a real efficiency improvement.
- **In areas where I have partial understanding**, AI can greatly shorten the time from 0 to 1 breakthrough. This would have been hard to imagine in the past, but this year, getting started with many new topics actually took only half a day.
- **In completely unfamiliar domains**, without sufficient background knowledge to verify and filter AI's output, one can only passively try each suggestion AI gives, which easily leads to spending a lot of time on fruitless efforts.

Simply put, my own ceiling determines AI's ceiling. In the end, the user still needs to study hard and improve themselves.
