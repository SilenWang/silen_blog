---
title: What AI Can and Cannot Do - Two Failed Attempts
categories: Others
date: 2025-07-15 00:44:45
tags: ['AI', 'aider', 'openfyde', 'theia']
---

Compared to two years ago at my previous company, I now use AI much more frequently. If the limitation back then was the convenience of the tools themselves, the current limitation seems to be my own knowledge and technical skills rather than AI.

<!-- more -->

Last weekend I made two attempts. The first was trying to use Aider and Openhands to help me modify Theia's remote module. According to AI's explanation, Theia's design is completely frontend-backend separated - the frontend only displays while the backend handles all requests. If this is correct, developing a remote version that works in the browser version of Theia should be feasible. With this, I could permanently solve my pain point when using VSCode on Fydetab Duo - insufficient graphics performance causing VSCode to lag.

However, although both agents could provide some guidance, without basic knowledge of Node.js/Electron, I couldn't determine whether the code written by the agents was reliable just based on their responses. Since I also didn't know how to debug a large Node.js project like Theia, I couldn't identify the issues from error messages. In the end, with my current knowledge and skills, I still couldn't complete the development under AI's guidance.

The second attempt was to solve the same problem but use a different method. Since I couldn't build my desired IDE, could I improve the graphics performance of Fydetab Duo? Specifically, by updating the kernel to enable Vulkan support on FydeOS running on Duo. With Vulkan, VSCode's performance might improve significantly.

This attempt also failed. While AI guided me to find the new kernel source code that might work with RK3588, I couldn't fully understand the descriptions about modifying/replacing content during compilation. So I tried the most brute-force approach - directly replacing the kernel project repo in overlay, forcibly using 6.12 as 6.1. I thought at worst some features wouldn't work? I underestimated it - it bricked my engineering version Duo...

Again, the same problem - with my knowledge and skills, I didn't know whether what I was doing was correct, and it was hard to get effective improvement information from the results. I only knew it didn't work, but how to make it work? Based on two years of experience, blindly trying AI's suggestions one by one isn't reliable, as there might be many solutions that professionals would immediately recognize as wrong. Without any foundation, after endless trials, I couldn't get any positive feedback.

Therefore, I think the AI I can currently use is truly just a Copilot. I shouldn't expect it to implement anything from scratch that I have no concept of. Only after establishing basic knowledge and being able to get effective feedback (at least being able to test and understand error messages) can I maximize its value as a Copilot.

That is, going from 0 to 1 is still difficult, but after getting started, going from 1 to 10 should be easier than before. As for 10 to 100... I've never reached 100, so I really don't know... But I think... these might be limited by context length.
