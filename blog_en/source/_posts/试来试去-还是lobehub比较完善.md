---
title: After Trying Many Options, LobeHub Still Stands Out
categories: Coding
date: 2026-04-06 00:32:41
tags:
  - AI Agent
  - LobeHub
  - MCP
  - Development Tools
---

Recently, I spent several days deeply experiencing multiple AI Agent products, from WorkBoddy, OpenClaw, LobeHub, DeerFlow2, LirbeChat, DingTalk Wukong to OpenCode, trying almost every tool I could find. Everyone has different needs, but my main focus is on **open source and self‑deployable**, **feature completeness**, and **practical usability**. Below is my personal comparison from a practical usage perspective.

<!-- more -->

## Feature Comparison

LobeHub (formerly Lobechat) has undergone its second major refactoring, so its features are quite rich:

| Feature | LobeHub | OpenClaw | DeerFlow2 | WorkBoddy | LirbeChat | DingTalk Wukong | OpenCode |
|---------|:-------:|:--------:|:---------:|:---------:|:---------:|:---------------:|:--------:|
| Open Source | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ | ✅ |
| Self‑deployable | ✅ | ✅ | ✅ | ❌ | ❌ | ❌ | ✅ |
| Complete Web Features | ✅ | Fair | ✅ | ❌ | ❌ | ❌ | ✅ |
| Multiple Search Support | ✅ | ✅ | ✅ | Limited | ✅ | Depends on DingTalk | ✅ |
| Third‑party LLM Provider | ✅ | ✅ | Limited | ❌ | ❌ | ❌ | ✅ |
| MCP Support | ✅ | ✅ | ❌ | ❌ | ❌ | ❌ | ✅ |
| Runtime Support | Cloud Sandbox | subprocess | Container | Local | Local | Local | Local |
| Scheduled Tasks | ✅ | ❌ | ❌ | ✅ (Enterprise) | ❌ | ✅ | ❌ |
| Agent Generation Assistant | ✅ | ✅ | ✅ | Built‑in Marketplace | Built‑in Marketplace | Built‑in Marketplace | ❌ |

## Current Limitations of LobeHub

- **No self‑hosted runtime**: The official team currently has no plans to develop a self‑hosted sandbox. Running through the cloud can be somewhat unstable, but theoretically this could be replaced via skills and MCP.
- **Scheduled tasks unavailable in self‑hosted version**: Scheduled tasks are disabled in the self‑hosted version. The official explanation is that this feature is provided via Qstash. In theory, one could write their own skills and MCP services, but there would be no corresponding monitoring dashboard.

## Final Thoughts

Objectively speaking, each product has its own suitable scenarios. OpenClaw is more developer‑friendly, DingTalk Wukong is more mature in enterprise settings, and OpenCode focuses on coding tasks. I needed a balance of **open source + self‑deployment + comprehensive features**, and LobeHub is currently the most suitable choice. However, it is not perfect—scheduled tasks and local runtime in self‑hosted scenarios remain weak points. I look forward to further updates in future self‑hosted versions.
