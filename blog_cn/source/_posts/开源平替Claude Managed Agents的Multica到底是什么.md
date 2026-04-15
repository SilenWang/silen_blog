---
title: 开源平替Claude Managed Agents的Multica到底是什么
date: 2026-04-15 16:00:00
tags:
- AI
- Claude
- 开源
categories: Coding / Others
---

前几天在量子位看到一篇报道，说Anthropic刚发布Claude Managed Agents就被开源项目"平替"了,这个开源项目就是Multica。刚好最近也在关注AI Agent的各种进展,就深入了解了一下,发现一些报道可能略有偏颇,这里做个记录。

<!-- more -->

## 时间线梳理

首先,报道说Multica是"被开源秒跟",这个表述并不完全准确。根据项目仓库的信息:

- Multica**今年2月**就正式创立了
- Claude Managed Agents是**今年4月**才宣布的

所以准确地说,Multica反而是更早出发的那个。4月的star暴涨确实是在Claude Managed Agents发布之后,但这更像是"英雄所见略同",而非单纯的跟随。

目前Multica已经来到了**v0.2.0**,本周刚刚更新,增加了自动化任务(cron)模块,开发速度相当惊人。

## Multica核心是什么

根据官方README,Multica给自己的定位是**开源的托管Agent平台**,核心理念是"把coding agents变成真正的队友"。

几个核心功能:

1. **Agent即队友** - 可以像给人分配任务一样给Agent分配issue,Agent会自动领取、写代码、报阻塞、更新状态
2. **全流程自主执行** - 配置好后免运维,支持任务排队、认领、执行、完结/失败的全生命周期
3. **Skills沉淀** - 每一套解决方案都会转化为可复用的Skill,团队共享
4. **统一算力运行时** - 一个控制台管控所有算力资源,兼容本地和云端
5. **多工作区隔离** - 按团队组织工作,工作区级别隔离

支持Claude Code、Codex、OpenClaw和OpenCode。

## 一些有趣的观察

1. 项目本身也是用Claude Code辅助开发的,不知道Anthropic作何感想
2. 创始人之前还创办了面向开发者的AI搜索引擎Devv.ai
3. 目前13.2k stars,1.6k forks,在开源项目中增长相当快
4. 技术栈是Next.js + Go + PostgreSQL,挺标准的现代架构

## 怎么看待

我觉得与其说Multica是"平替",不如说它代表了**一种需求趋势** - 大家都希望能更系统地使用AI Agent,而不仅仅是把它当作一个孤立的对话工具。Claude Managed Agents解决的是企业级托管的问题,Multica解决的是开源可自部署的问题,两者定位有所差异。

如果你也对多Agent协作、团队级AI赋能有兴趣,可以关注一下这个项目。官方给出的安装方式很简洁:

```bash
# macOS/Linux
brew install multica-ai/tap/multica

# 或者
curl -fsSL https://raw.githubusercontent.com/multica-ai/multica/main/scripts/install.sh | bash

# 然后一行命令完成配置
multica setup
```

感兴趣的朋友可以自行体验。后续如果有什么使用心得,再分享。