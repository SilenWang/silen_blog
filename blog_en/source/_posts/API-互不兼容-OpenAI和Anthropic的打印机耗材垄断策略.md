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

Lately I've been working on integrating several AI applications, which means I need to switch between different models. That's when I discovered a truly infuriating fact – the APIs of OpenAI and Anthropic are **almost completely incompatible**.

<!-- more -->

## Starting with a Simple Model Switch

Suppose you're already using OpenAI's GPT‑4o and now you want to try Claude 3.5 Sonnet. Intuitively, you'd just change the endpoint and the model name, right? Naïve.

```python
# OpenAI calling style
import openai
client = openai.OpenAI(api_key="sk-xxx")
response = client.chat.completions.create(
    model="gpt-4o",
    messages=[{"role": "user", "content": "Hello"}]
)
```

And Anthropic?

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

At first glance they look similar, but once you actually use them, pitfalls are everywhere:

- **Different client SDKs**: one uses the `openai` package, the other uses `anthropic` – you can't mix them
- **Different parameter names**: `max_tokens` appears in different places and with different syntax
- **Different response formats**: OpenAI streams with `for chunk in response`, Anthropic uses `for event in response` and you have to check the event type
- **Different streaming implementation details**: even handling streaming responses has to be rewritten
- **Different tool/function calling protocols**: this is the killer – if your application uses function calling, switching models essentially means rewriting the entire invocation layer
- **Different system prompt formats**: OpenAI uses a `system` role, Anthropic has its own system prompt parameter
- **Different image input formats**: even sending an image requires a different format

These differences force developers either to write a ton of adapter code or to use third‑party libraries (like LiteLLM, OpenRouter) as a wrapper.

## This Reminds Me of Printers

The situation is all too familiar. Back when you bought a printer, the machine itself could be sold cheaply, but the ink/toner cartridges were always proprietary. HP printers only accept HP cartridges, Canon only Canon. Third‑party compatible cartridges are either blocked by firmware or stopped by "official authentication chips".

Think about it this way:
- **OpenAI = HP**: big, well‑established ecosystem, everyone uses its API format
- **Anthropic / Claude = Canon**: a later entrant with good performance, but insists on its own protocol
- **You = the locked‑in user**: once you deeply integrate with one API, the migration cost becomes terrifying

## The Three Traps of API Lock‑In

The incompatibility between OpenAI and Anthropic isn't purely "technical debt" or "different design philosophies" – there's a clear business logic behind it:

1. **Switching cost**: once you've written a lot of code and deeply integrated with one API, it's hard to just switch. Even if another model is better or cheaper, the psychological and engineering migration costs remain.
2. **Ecosystem lock‑in**: a huge ecosystem has grown around the OpenAI API – various tools, monitoring, proxies, caching layers. Switching APIs means all of these have to be adapted too.
3. **Data moats**: although the API doesn't directly lock your data, when you deeply rely on a platform's specific features (e.g., OpenAI's structured output, Anthropic's extended thinking), you become even more tied down.

Isn't this the digital‑age version of printer ink? The machine (API standard) itself can be open, but the consumables (implementation details, ecosystem tools, proprietary features) are proprietary.

## Fortunately, There Are Some Disruptors

Of course, the community isn't sitting idle. Quite a few projects are trying to break this lock‑in:

- **LiteLLM**: provides a unified interface covering 100+ model providers
- **OpenRouter**: one API to call multiple models
- **Vercel AI SDK**: a framework‑level abstraction layer
- **LangChain**: although often criticised for being over‑engineered, it's genuinely useful in this scenario

But essentially all of these adapt to the OpenAI API style – Anthropic's native API hasn't become a "standard" because of them.

## There Is Nothing New Under the Sun

From printers to smartphone charging ports, from game‑console exclusives to cloud‑provider lock‑in strategies, commercial companies have been using incompatibility to hold users hostage for over a hundred years.

The AI era is no exception. The incompatibility between OpenAI and Anthropic's APIs is, plain and simple – **the old printer‑vendor trick of "buy the machine freely, but you must use my consumables" has just put on a new coat**.

As developers, what we can do is:
1. Reserve an abstraction layer when designing systems – don't tie your code to any single API
2. Keep an eye on standardisation progress, e.g., whether a unified API standard emerges
3. Support projects that are committed to an open ecosystem

But honestly, expecting commercial companies to voluntarily give up lock‑in strategies is like expecting printer manufacturers to give away ink cartridges for free – dream on.

There is nothing new under the sun.
