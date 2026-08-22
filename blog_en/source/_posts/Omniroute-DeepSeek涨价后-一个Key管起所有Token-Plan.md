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

DeepSeek raised its prices again recently. As someone whose token usage is closing in on 100M tokens a month (see {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [My Token Usage Exploded 13× in Three Months] %}), that news hit me hard: I can't cut down usage, so I had to get clever about the unit price.

<!-- more -->

## DeepSeek's Price Hike, and a 100M+ Token Monthly Bill

Some background first. When DeepSeek V4-Pro launched, it announced a permanent price cut to a quarter of the original price, and I even wrote a post praising them for it ({% post_link DeepSeek-永久降价-大模型价格战的新玩法 [DeepSeek Finally Changed Its Pricing... Wait... A Price Cut?] %}). A few months later, the prices went back up...

And my usage is no longer anywhere near that original scale. In that earlier post I calculated that my token consumption exploded from 38M to 500M per month in three months, and lately it has settled in at over 100M a month. Usage went up, unit price went up — the math just doesn't work anymore.

## The "Sweet Spot" Tier Paradox of Token Plans

To save money, I started looking into vendors selling Token Plans (monthly subscription / prepaid token packages), and I noticed a very interesting pattern:

> Almost every vendor has a "sweet spot" tier — by per-token unit price, it's usually the cheapest one; the more expensive packages actually cost more per token.

This is probably deliberate: the cheap tier is a loss leader to attract customers, while the large packages are aimed at users who "can't be bothered," so they carry a higher unit price. In other words, **blindly buying the most expensive package is actually the worst deal**.

So the right way to save money is clear:

1. Pick multiple vendors and buy each one's "sweet spot" tier
2. Register multiple accounts to add up enough total quota
3. After all, models from DeepSeek, Kimi, GLM, and MiniMax are roughly comparable — use whichever is cheapest

But this approach has one huge problem: **too many keys**. Say you have 5 vendors × 2 accounts = 10 keys, and your toolchain includes Claude Code, Codex, opencode... Every tool needs its key configured, keys swapped, and you have to remember which vendor has how much quota left. Managing all those keys alone is enough to drive you crazy.

## Omniroute: One Key to Rule Them All

That's where a project like Omniroute comes in.

[Omniroute](https://github.com/diegosouzapw/OmniRoute) is an open-source AI gateway (MIT license, with over 50k stars on GitHub). What it does in one sentence: **it aggregates hundreds of model providers into a single unified endpoint, so the client only needs one key**.

Key features:

- **Huge provider coverage**: 340+ providers (90+ of which offer free tiers) and 1200+ models, covering all the mainstream ones — Kimi, Claude, GPT, Gemini, GLM, DeepSeek, MiniMax
- **Rich aggregation strategies**: 19 built-in routing strategies, four-tier auto-fallback (Subscription → API → Cheap → Free), quota-aware auto-fallback, circuit breakers, and key cooldowns
- **Zero changes on the client side**: just point the `base_url` of Claude Code / Codex / Cursor to `http://localhost:20128/v1`; the tools only see one key, and everything else is switched in the Omniroute dashboard
- **It even saves tokens**: a built-in RTK + Caveman compression engine that claims to save 15%–95% of tokens
- **Flexible deployment**: npm global install, Docker, an Electron desktop app, or PWA

The general flow: deploy an Omniroute instance (I used Docker), fill in each vendor's key in the dashboard, configure routing strategies (e.g. "prefer the cheapest one, auto-fallback when quota runs out"), then point all your AI tools' `base_url` at it. Want to switch vendors afterward? One click in the dashboard, completely transparent to the client.

Omniroute has also made a lot of usability improvements. For example, it ships an **MCP server** (exposing over 100 tools), so AI can configure things for you — switching combos, setting routing strategies, checking quota, all by just asking. For a lazy person like me, that saves a lot of hassle.

## Pitfalls I Hit in Real Deployment

That said, Omniroute is far from "works out of the box" in practice. I hit two pretty annoying pitfalls.

### Pitfall 1: The search interception feature is buggy — Codex and Claude Code can't search the web

Omniroute has a "search interception" feature: it intercepts the client's native search tool (e.g. Claude Code's `web_search_20250305`), forwards it to its own `/v1/search` gateway, and runs it against the search provider you configured (Exa, Serper, etc.).

Great idea, but there are bugs. Claude Code sends the versioned `web_search_20250305`, while Omniroute's interception matcher only matches the bare `web_search` — the version suffix doesn't match, so interception never kicks in, the upstream model returns empty search results, and Claude Code ends up "searching into the void." Codex has a similar problem: its standalone search endpoint isn't supported for custom providers.

This pitfall hurts a lot: the whole point of aggregation is "one key for everything," but the search feature just dies, effectively cutting the AI tools' web access in half. The relevant issues are still open ([#9772](https://github.com/diegosouzapw/OmniRoute/issues/9772), [#9725](https://github.com/diegosouzapw/OmniRoute/issues/9725), [#8674](https://github.com/diegosouzapw/OmniRoute/issues/8674)), so for now we just wait for fixes.

### Pitfall 2: Very new providers can't sync quota

Omniroute's quota sync works on an allowlist — each provider's quota-checking endpoint has to be registered individually, and new providers frequently get missed.

For example, the Qwen Cloud Token Plan I use works fine for inference, but its quota never shows up in the dashboard's Provider Quota. Worse, when quota can't be synced, the "quota-aware auto-fallback" breaks: exhausted accounts aren't marked correctly, the multi-account rotation gets messed up, and in the worst case every account gets labeled "quota exhausted" and never recovers automatically after the reset ([#9603](https://github.com/diegosouzapw/OmniRoute/issues/9603)).

In other words, **the newer and cheaper the provider, the more likely you'll hit quota sync problems** — and those are exactly the ones we most want to use.

## Summary

After using it for a while, here's my take:

- Omniroute solves the core problem — **one key for all accounts** — and the aggregation strategies genuinely work well and save a lot of hassle
- But it's still young: **peripheral features like search interception and quota sync have quite a few pitfalls**, and newer providers are more likely to hit them
- If you just use it as a "multi-key aggregator," the experience is great; but if you rely on its search gateway or expect every provider to sync quota, check the status of the relevant issues first

As usual, here's the project link: [github.com/diegosouzapw/OmniRoute](https://github.com/diegosouzapw/OmniRoute). Feel free to dig in if you need it.
