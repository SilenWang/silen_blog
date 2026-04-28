---
title: DeepSeek-V4 is Finally Here
categories: AI
date: 2026-04-28 13:00:00
tags:
  - DeepSeek
  - AI
  - Large Language Model
  - Reasoning
  - Agent
---

A few days ago, DeepSeek officially released the V4 preview, and this version really attracted a lot of attention. Rumors had been circulating since before the Chinese New Year. Since I happened to be leaving my job at the time, I didn't get to try it out immediately. Now that I have some time, I integrated it with Opencode and Multica and wrote this blog post to test it out.

<!-- more -->

## API Availability

So far, the DeepSeek API is accessible at any time of day without experiencing peak-hour queuing issues during specific time periods. They launched an unprecedented price reduction plan (who else slashes prices by 90%... I'm impressed), which also suggests that this release hasn't caused the kind of unprecedented frenzy as before, and their computing power probably isn't fully utilized.

## Agent Capability Optimization

This is probably the most interesting part of V4 for me. DeepSeek has specifically adapted and optimized for mainstream Agent products including **Claude Code, OpenClaw, OpenCode, CodeBuddy**, and others. In terms of coding tasks and document generation, while I wouldn't say the effect has improved dramatically, the stability has certainly increased significantly.

The previous R1, although quite good at scoring and daily responses, always had some inexplicable issues with tool calling. That's why in both Aider and Opencode, I used the non-thinking V3 to ensure normal tool calling and prevent the Agent from getting stuck for no reason.

## Consumption Testing

Previously, when using aider-chat + DeepSeek, I translated dozens of articles in my blog project, using bash scripts for task control. The entire translation took about an hour. With this DeepSeek update, I used opencode + DeepSeek, letting the Agent itself control the tasks directly, to see how well it worked and how long the execution took.

In actual testing, a simple prompt (`Clone silen_blog, create a new branch, and in that branch, gradually review all Chinese blog md files in the posts folder under the blog_cn directory, check for inappropriate Chinese expressions and typos, make corrections, and after completing the corrections, push the changes back to GitHub`) let the Agent proofread my articles for Chinese expression in just 30 minutes, with commits affecting 105 files, while I have about 150+ MD blog files. Since it was just a simple test, I didn't verify each blog one by one to see if they were all actually reviewed and modified, but from the Diff results, it did find quite a few typos.

In terms of price, using Flash cost only 1.5 RMB, which is quite cheap—especially considering this was through Opencode, where token consumption theoretically grows exponentially compared to Aider's previous usage. However, for full-text translation, the consumption is only a few yuan.

## Personal Expectations

Currently, many of my work and ideas depend heavily on AI. So if I want to work but can't use AI because of queuing, it's really uncomfortable... However, amidst the AI frenzy in China, even Minimax, an AI with fewer parameters, experiences stable request blockages at 3 PM... I hope that in the second half of the year, DeepSeek will not only reduce API prices but also make it possible for individuals to run a reasonably functional AI at a relatively low cost.
