---
title: AI Agent Makes Open Source Software Easier to Leverage
categories: AI
date: 2026-06-03 22:00:00
tags:
  - AI Agent
  - Open Source
  - BananaSlice
  - Programming
  - LLM
---

Recently I submitted a PR to BananaSlice ([#11](https://github.com/IrfanulM/BananaSlice/pull/11)), adding configurable Gemini-compatible API endpoint support. The whole process made me realize: AI Agents not only help people write code faster, but also change how we leverage open source software.

![banner](https://raw.githubusercontent.com/silenwang/Gallary/master/2026/06/upgit_20260615_1781454783.webp)

<!-- more -->

## The Old Problem: So Close, Yet So Far

I think many people have experienced this: you find an open source tool with an appealing feature set, but it's always missing one or two small features to fully fit your use case. Take BananaSlice for example — it's a tool that uses Gemini's image API for quick image editing, but it only supports the official Gemini API endpoint. If you want to use your own proxy or a third-party compatible service (like aihubmix), you're out of luck.

In the past, when faced with this situation, you'd either fork the project and dig into the code yourself, or switch to a different tool. The problem is — open source projects use all kinds of programming languages, and no one can be proficient in all of them. Learning a new language or diving deep into an unfamiliar codebase just for a "one little thing" is simply not worth the investment.

## What AI Agent Changed

The emergence of AI Agents has fundamentally changed this dilemma. When the required changes are small, you can just use AI to implement them without needing to deeply understand that unfamiliar language or framework.

Take BananaSlice again: its core tech stack is Rust + TypeScript. I can read some TypeScript, but I'm not a Rust developer. If this were two years ago, I'd probably have given up — spending time learning Rust just for a simple config feature is far more costly than the value of the feature itself.

But this time I used Opencode + DeepSeek, described the requirements (support for custom Base URL and Model ID), and the Agent automatically handled code modification, frontend UI addition, data persistence, and everything else. All I had to do was review whether the generated code was relevant to my needs and test if the feature worked properly.

From submission to merge by the project owner, the entire process went smoothly. The project owner replied: "Always wanted to add this but never had the time to." — This is the norm for open source projects: maintainers have limited bandwidth, and many good improvements get shelved because there's no one to implement them.

## More Than "Making It Work" — It's Contribution

What's even more interesting is that AI Agents merge "fixing something for yourself" with "contributing to the project."

Once the changes are done, it's natural to submit a PR back to the original project. For open source projects, this means:
- **Users** get the features they need quickly
- **The project itself** improves and becomes more feature-complete
- **All users** benefit from the enhancement

Those "almost there" gaps that were previously abandoned due to language barriers or learning costs can now be filled quickly with AI, while giving back to the community and helping software iterate faster.

## Summary

AI Agents have lowered the barrier to leveraging open source software to an unprecedentedly low level. You no longer need to be proficient in a language to participate in or modify a project. As long as the requirements are clear and the changes are reasonable, AI can help you through the entire process from code implementation to PR submission.

This not only improves individual productivity but also accelerates the evolution of the entire open source ecosystem. Next time you find an open source project you like but is missing one or two features, give AI Agent a try — maybe your PR is exactly the push the project needs.
