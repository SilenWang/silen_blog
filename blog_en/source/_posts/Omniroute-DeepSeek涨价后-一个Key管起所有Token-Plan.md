---
title: "Omniroute: After DeepSeek Raised Prices, I Managed All My Token Plans with One Key"
categories: AI
date: 2026-08-22 22:00:00
tags:
  - Omniroute
  - DeepSeek
  - Token
  - AI
  - API
  - Cost-saving
---

DeepSeek recently officially raised its prices. As someone whose token usage was closing in on 100M tokens last quarter (see {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [My Token Usage Exploded 13× in Three Months] %}), I didn't pay much attention at first — after all, DeepSeek was already very cheap, and even the raised price was still lower than the official price announced when V4 first launched. But I underestimated my usage: after a day of real-world testing at the new prices, it became clear that the cost is genuinely hard for an individual to bear. If I didn't find a way around it, I'd be working just to pay for tokens.

<!-- more -->

## The "Sweet Spot" Tier Paradox of Token Plans

To save money, I started looking into vendors selling Token Plans (monthly subscription / prepaid token packages), and I noticed a very interesting pattern:

> Almost every vendor has a tier that offers the best value — by per-token unit price, it's the most cost-effective; while the more expensive packages actually cost slightly more per token.

This is probably deliberate: the cheap tier is a loss leader to attract customers, while the large packages are aimed at users who "can't be bothered," so they carry a higher unit price. In other words, **blindly buying the most expensive package isn't necessarily a good deal**.

So the right way to save money is clear:

1. Pick multiple vendors and buy each one's best-value tier
2. Register multiple accounts to add up enough total quota
3. After all, DeepSeek, Kimi, GLM, and MiniMax all offer these models — use whichever is cheapest

But this approach has one huge problem: **too many keys**. Say you have 5 vendors × 2 accounts = 10 keys, and your toolchain includes Claude Code, Codex, opencode... Every tool needs its key configured, keys swapped, and you have to remember which vendor has how much quota left. Managing all those keys alone is enough to drive you crazy.

## Omniroute: One Key to Rule Them All

That's where a project like Omniroute comes in.

[Omniroute](https://github.com/diegosouzapw/OmniRoute) is an open-source AI gateway (MIT license, with over 50k stars on GitHub). What it does in one sentence: **it aggregates hundreds of model providers into a single unified endpoint, so the client only needs one key**.

It has detailed feature comparisons on its official site; in short, it offers broad provider support and rich routing strategies.

The general flow: deploy an Omniroute instance (I used Docker), fill in each vendor's key in the dashboard, then use an agent to call the MCP server to configure the routing strategies (e.g. "prefer the cheapest one, auto-fallback when quota runs out"), and point all your AI tools' `base_url` at it. Want to switch vendors afterward? One click in the dashboard, completely transparent to the client.

Omniroute has also made a lot of usability improvements. For example, it ships an **MCP server** (exposing over 100 tools), so AI can configure things for you — switching combos, setting routing strategies, checking quota, all by just asking. For a lazy person like me, that saves a lot of hassle.

## Pitfalls I Hit in Real Deployment

That said, Omniroute isn't entirely problem-free in practice.

### Pitfall 1: The search interception feature is buggy — Codex and Claude Code can't search the web

Omniroute has a "search interception" feature: it intercepts the client's native search tool (e.g. Claude Code's `web_search_20250305`), forwards it to its own `/v1/search` gateway, and runs it against the search provider you configured (Exa, Serper, etc.).

The feature itself is actually great — both Codex and Claude Code need server-side search results, and many overseas providers only sell tokens, not tools like that. But it can't actually be used yet. Claude Code sends the versioned `web_search_20250305`, while Omniroute's interception matcher only matches the bare `web_search` — the version suffix doesn't match, so interception never kicks in, the upstream model returns empty search results, and Claude Code ends up "searching into the void." Codex has a similar problem: its standalone search endpoint isn't supported for custom providers.

This pitfall's impact is neither huge nor trivial — not huge, because you can get results another way, e.g. letting the agent search via skills or MCP; but not trivial either, because the plan was for Omniroute to solve everything, and now we just have to wait for fixes. The related issues are still open ([#9772](https://github.com/diegosouzapw/OmniRoute/issues/9772), [#9725](https://github.com/diegosouzapw/OmniRoute/issues/9725), [#8674](https://github.com/diegosouzapw/OmniRoute/issues/8674)).

### Pitfall 2: Very new providers can't sync quota

For example, the Commandcode Go I use works fine for inference, but its quota never shows up in the dashboard's Provider Quota. Worse, when quota can't be synced, the "quota-aware auto-fallback" may break: exhausted accounts aren't marked correctly, the multi-account rotation gets messed up, and in the worst case every account gets labeled "quota exhausted" and never recovers automatically after the reset ([#9603](https://github.com/diegosouzapw/OmniRoute/issues/9603)).

In other words, **the newer and cheaper the provider, the more likely you'll hit quota sync problems** — and those are exactly the ones we most want to use.

### Pitfall 3: Concurrency is capped by default — your agent suddenly drops mid-task

The first two pitfalls are tolerable, but this one is genuinely infuriating: **out of the box, Omniroute doesn't allow heavy concurrent tasks, and your agent will abruptly die mid-run.**

Modern agents (Claude Code, Codex, opencode, etc.) never fire one request at a time. A big task spawns several subtasks in parallel, each blasting out a batch of LLM requests — subagents in opencode, parallel tool calls in Claude Code. Concurrency spikes instantly.

On the Omniroute side, **every account (connection) has a concurrency cap by default**, and requests beyond it just queue up: the account semaphore queues up to 20 by default with a 30-second wait, then starts rejecting; the request queue has its own execution limits too. When a heavy task pushes concurrency past the cap, the overflow requests either sit queued or get rejected outright. Many client agents treat these errors (429, timeouts, queue-full) as fatal and kill the whole session — which is exactly what you see as "suddenly dropped mid-use."

What's worse, the limit applies **per account (connection), not globally**. Think you've set a global concurrency cap? No — it only constrains that one account; other accounts and other providers each have their own quota. For a true global cap, you'll just have to wait ([#7778](https://github.com/diegosouzapw/OmniRoute/issues/7778)).

My workarounds:

1. Set `max_concurrent` on each connection to match the upstream account's real concurrency (e.g. 1–2 for providers like GLM that can only handle one or two concurrent requests). Don't leave them unbounded by default, and don't get greedy with high values.
2. When you really need to run heavy tasks in parallel, either cap the concurrency on the agent side, or bump the queue settings under `Settings → Resilience → Request Queue`.
3. When debugging, watch out for session affinity pinning your session to a saturated account and making it idle-wait the full 30 seconds.

## Wrapping Up

It's not perfect, but this setup covers at least 80% of my needs for now. Hopefully those issues get fixed soon...

As usual, here's the project link: [github.com/diegosouzapw/OmniRoute](https://github.com/diegosouzapw/OmniRoute). Feel free to dig in if you need it.
