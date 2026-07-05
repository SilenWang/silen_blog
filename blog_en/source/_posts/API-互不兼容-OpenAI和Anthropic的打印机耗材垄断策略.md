---
title: API Incompatibility – OpenAI and Anthropic's "Printer Ink Monopoly" Strategy
date: 2026-07-01 22:00:00
tags:
- AI
- API
- OpenAI
- Anthropic
- Rant
categories: Others
---

Lately I've been working on integrating several AI applications, which means I need to switch between different models. That's when I discovered a problem — the APIs of OpenAI and Anthropic are **incompatible**.

![banner](https://raw.githubusercontent.com/silenwang/Gallary/master/2026/07/upgit_20260706_1783268947.webp)
<!-- more -->

## Starting with a Simple Model Switch

Suppose you're already using OpenAI's GPT‑4o and now you want to try Claude 3.5 Sonnet in the same tool. Normally, you'd just change the endpoint and the model name, right? Not really. Here's the official OpenAI SDK call:

```python
# OpenAI calling style
import openai
client = openai.OpenAI(api_key="sk-xxx")
response = client.chat.completions.create(
    model="gpt-4o",
    messages=[{"role": "user", "content": "Hello"}]
)
```

This was the protocol that everyone basically followed after ChatGPT took off. But as the technology evolved rapidly, Anthropic introduced their own way of doing things:

```python
# Anthropic calling style
import anthropic
client = anthropic.Anthropic(api_key="sk-ant-xxx")
response = client.messages.create(
    model="claude-3-5-sonnet-20241022",
    max_tokens=1024,
    messages=[{"role": "user", "content": "Hello"}]
)
```

At first glance they look similar, but actually dealing with both is quite troublesome:

- **Different client SDKs**: one uses the `openai` package, the other uses `anthropic` – you can't mix them
- **Different parameter names**: `max_tokens` appears in different places and with different syntax
- **Different response formats**: OpenAI streams with `for chunk in response`, Anthropic uses `for event in response` and you have to check the event type
- **Different streaming implementation details**: even handling streaming responses has to be rewritten
- **Different tool/function calling protocols**: this is the killer – if your application uses function calling, switching models essentially means rewriting the entire invocation layer
- **Different system prompt formats**: OpenAI uses a `system` role, Anthropic has its own system prompt parameter
- **Different image input formats**: even sending an image requires a different format

I guess this is exactly why DeepSeek created an Anthropic‑compatible endpoint...

## This Reminds Me of Printers

This scenario is all too familiar. Back when you bought a printer, the machine itself could be sold cheaply, but the ink/toner cartridges were always proprietary. Brand A's printer only works with Brand A's cartridges, Brand B's only with Brand B's. Third‑party compatible cartridges are either blocked by firmware updates or stopped by "official authentication chips".

I think the incompatibility between OpenAI and Anthropic's APIs isn't purely "technical debt" or "different design philosophies" – there's a clear business logic behind it:

1. **Switching cost**: once you've written a lot of code and deeply integrated with one API, it's hard to just switch. Even if another model is better or cheaper, the psychological and engineering migration costs remain.
2. **Data moat**: although the API doesn't directly lock your data, when you deeply rely on a platform's specific features (e.g., OpenAI's structured output, Anthropic's extended thinking), you become even more tied down.

Isn't this the digital‑age version of printer consumables? The machine (API standard) itself can be open, but the consumables (implementation details, ecosystem tools, proprietary features) are proprietary.

## There Is Nothing New Under the Sun

From printers to smartphone charging ports, from game‑console exclusives to cloud‑provider lock‑in strategies, commercial companies have been using incompatibility to hold users hostage – this trick is nothing new.

As a developer and a user, besides complaining... I still have to use them anyway. Sigh.
