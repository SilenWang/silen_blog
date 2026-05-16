---
title: "AI Tools Experience & Rant: Multica, MiniMax, Gemini, Notebooklm"
date: 2026-05-16 07:30:00
tags:
  - AI
  - Multica
  - MiniMax
  - Google One
  - Notebooklm
  - Gemini
  - rant
categories: AI
---

I recently switched to another batch of commonly used AI tools — some pleasantly surprising, some left me speechless. Here's my honest record.

## Multica: Master of Cyber Livestock Management

Multica's recent updates are genuinely impressive. The development pace has always been fast. Compared to last time I mentioned it, they've now added **scheduled tasks** and **squad features** — you can directly set timed triggers to let agents work, and collaborate with team members.

For indie developers, Multica's approach as an "AI-native task management platform" is indeed ahead of its time. With just one-time configuration, you no longer need to switch back and forth between multiple tools — you can set up, assign, and execute tasks through a unified platform. At the same time, lookback/review is also more convenient than directly browsing Agent software records.

## MiniMax: Gradually Showing Its Disadvantage

At my previous company, I always used MiniMax's Plan tier — cheap, enough quota, and could complete the tasks I needed. But recently I've gradually felt it's not as good, especially after DeepSeek V4 came out — the contrast is obvious.

First, **context length**. To be honest, the scripts I deal with are all relatively short. In theory, I don't need very long context. But when I want the Agent to fully autonomously perform a round of debugging, because it needs to read a lot of runtime logs, MiniMax frequently exceeds the length limit and exits... With DeepSeek V4's million-token context, I can't say it gets it right on the first try, but at least it can complete the task, and then I can investigate remaining issues based on the execution logs.

Second, **answer accuracy**. This is admittedly quite subjective, but in my recent deep learning environment configuration process, the probability of MiniMax getting it right on the first try is indeed not high. DeepSeek, on the other hand — although the environment it configures may not always be suitable — at least it can run.

Of course, models can improve over time, so I'll keep an eye on it. But in the near future, I probably won't continue using MiniMax much...

## Google One: Surprisingly Worth It

This one I actually want to praise. The Google One subscription now bundles **NotebookLM** and **multimodal Gemini** access. For anyone who reads academic papers regularly, it's a godsend.

NotebookLM handles various input types exceptionally well. Unlike some RAG tools I've used before, it can very accurately answer specific questions about a paper based on the uploaded information. At the same time, Google's internal multimodal integration allows me to draw schematic diagrams directly through NotebookLM, which helps me understand papers much faster (deep learning papers are full of formulas that are tough to parse). One subscription covering paper writing, literature review, and research — if you have a Japanese region account, the lowest tier is even cheaper than domestic vendors' memberships... The value proposition is solid.

## Gemini-cli: A KPI Product?

Finally, I have to rant about Gemini-cli. As mentioned earlier, Google's NotebookLM and Gemini models themselves are both excellent, but the CLI tool they made is just too perfunctory.

**Tool calling** is a poor experience — calls frequently go unanswered or return malformed output. **Output speed** is also slow; compared to similar streaming tools, it's noticeably slower. The whole experience feels like a mediocre teammate dragging down the team's star player. It makes you wonder whether this thing is a KPI product — "they have one, so we need one too." Not to say it's good — it's not even close to being usable.

---

These are my honest takeaways from recent AI tooling adventures. Some hits, some misses. Hope they all keep improving.
