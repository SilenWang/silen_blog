---
title: "My AI Infrastructure: After DeepSeek Raised Prices, I Built a Workspace That Doesn't Depend on Its API"
categories: AI
date: 2026-08-31 00:30:00
tags:
  - AI
  - DeepSeek
  - Gateway
  - Omniroute
  - New-API
  - APISIX
  - Bifrost
  - Codex
  - Claude Code
  - Self-hosted
---

This is the full story of how I rebuilt my AI workspace: why I decided I had to stop depending on DeepSeek's API, which options I looked at, and what problems each layer left me with. It's the sequel to {% post_link Omniroute-DeepSeek涨价后-一个Key管起所有Token-Plan [Omniroute: One Key to Rule All My Token Plans] %} — that post was about "I picked omniroute"; this one is about "after picking it, I found I had to stack several more layers on top".

<!-- more -->

## The price hike was just the trigger

When DeepSeek officially raised its prices, my first reaction was a shrug — it was already cheap, and even after the hike it was cheaper than the pricing announced when V4 first launched. Then I measured a single day at the new rates and got my math wrong: my volume is absurd (see {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [my token usage exploding 13× in three months] %}), so any percentage hurts, and if I didn't do something I'd be paying for the privilege of working.

But halfway through hunting for a cheaper option, I realised money wasn't the real problem.

My entire working environment — the container, the agents, the model config, the launch scripts — hung off **one cloud API endpoint**. A price hike I could live with. But what about the day they change the protocol, tighten rate limits, rename a model, or just have a bad afternoon? Most of the "my agent suddenly stopped working" incidents I've hit over the past couple of years came from that one layer, and I had zero control over it.

So the goal became: **build a working environment that doesn't depend on DeepSeek's API**.

Let me be precise, because this is easy to misread. I don't mean running models myself. I did look at local deployment seriously, but the machine I have can't run a model good enough to do real work, so the arithmetic killed it fast. What I wanted back was control of the **entry layer** — endpoint, routing, billing, protocol translation. Whoever actually serves the tokens can change; those four things cannot.

## What the thing looks like

Seen whole, my infrastructure is three layers, and in hindsight they can't be built in any other order:

```
container (my workspace, where the agents run)
  └─ launcher (a small Go program that decides where requests go this round)
      └─ gateway layer (omniroute / new-api / CLIProxyAPI / official direct, composable)
```

The bottom two existed long before this story. The gateway layer is what grew only recently. I even hand-wrote a tiny bridge program (`codex-deepseek-bridge`) early on, purely to translate between Codex and DeepSeek — that should have told me the protocol bill would come due eventually.

## Generation 1: straight to the official API

The first version was blunt: point the container config straight at DeepSeek's official API. It exposes two interfaces, and you use whichever your agent speaks — **Anthropic mode** for Claude Code and friends, **OpenAI mode** for Codex and friends (getting Codex onto DeepSeek at all was originally blocked by exactly this incompatibility; I ended up doing it the official way).

The config layer also solved one thing properly: a **vision fallback model**. When the main model can't see images, something has to catch that, so vision gets its own configurable fallback via environment variables.

Generation 1's problems weren't functional. They were all about config:

- **base_url scattered everywhere**. One copy for Codex, one for Claude Code, one for opencode, another for Reasonix — every agent added means one more place to edit.
- **Every new model has to be written into every one of those places.** When the gateway added a free model, I had to answer "did you update your config?" one client at a time.
- At some point I did a big flip: switch all traffic back to the official API, change the default model to a vision experimental build (keeping the stable-looking name in the UI, with the vision fallback pointed at the same one). The day I counted how many places that touched — that's when I knew this couldn't continue.

The one thing generation 1 got right was the Go launcher. Multiple forward targets, interactive pick at startup, configurable command, and it runs outside a container too. It was quietly already a router — it just knew how to switch everything at once, not how to split per request.

## Generation 2: omniroute, one key for everything

So I went to omniroute: aggregate every vendor's token plans behind a single entry point, one key in the client, switch vendors in a console while every client stays oblivious. I wrote up the whole thing in {% post_link Omniroute-DeepSeek涨价后-一个Key管起所有Token-Plan [the omniroute post] %}, so here are just the names of three pitfalls plus one new one.

1. Search interception matches bare tool names, not versioned ones — Codex and Claude Code can't search the web;
2. Brand-new vendors don't sync quota, so "quota-aware automatic fallback" silently dies;
3. Heavy concurrent tasks aren't allowed by default, so the agent just drops mid-session when sub-agents run in parallel.

The new one only surfaced later: **model metadata still has to be hand-copied**. A model existing in the gateway doesn't mean the clients can use it — Codex's `config.toml`, opencode's config, Reasonix's config each need their own edit. The gateway unified the key; it did not unify the model list.

More fundamentally: omniroute settled "which vendor do I send to" and never touched "which protocol do I send it with". And my agents don't speak one protocol.

## Generation 3: new-api, blocked by protocol

I moved to new-api expecting it to take "unify multiple protocols" further, which was exactly the disease.

It died on one very specific point: **upstream new-api doesn't support OpenAI's `responses` protocol**, and Codex speaks precisely that dialect. You configure the direct address locally, and the request is simply not supported. No config workaround — the two protocol stacks don't overlap.

## Generation 4: add CLIProxyAPI just to translate

The workaround was **another layer**. Install CLIProxyAPI to translate Codex's `responses` protocol and forward it to new-api.

```
Codex → CLIProxyAPI(:8317) → new-api → vendors
```

Its default port is 8317. On the Codex side you add a custom provider in `~/.codex/config.toml` with `base_url` pointed at it and `wire_api = "responses"` — that last one is Codex's protocol switch, and the entire reason this layer exists.

My chain at this point:

```
agents in the container (Codex / Claude Code / opencode / Reasonix / dsh / codebuddy)
   → CLIProxyAPI (protocol translation)
   → omniroute or new-api (routing and quota)
   → vendors
```

It works. And every failure became hard to place: when an agent misbehaves it could be the config in the container, CLIProxyAPI, the gateway, or the upstream vendor. Four layers, each with its own config and its own defaults. A large share of my debugging time went into establishing *which layer* was broken — and the method was cutting the chain shorter, one layer at a time.

## bifrost and apisix: evaluated, never deployed

Both were on my shortlist. I never actually deployed either, so what follows is "why I didn't ship it", not a pitfall log — I'd rather say that plainly than invent details.

**bifrost** pitches itself as a high-performance LLM gateway written in Go, and performance is exactly the right thing to care about when you burn this many tokens. The question is what it can do *at the protocol layer* for the mix I actually run — I need one entry that feeds both OpenAI-style and Anthropic-style clients. That's precisely where gateways differ most, and you can't tell from a README; you have to put it under your own traffic. At the time I evaluated it, omniroute was already carrying my daily load, so the evaluation stayed an evaluation.

**apisix** is a different animal entirely: a general-purpose API gateway with AI capabilities bolted on as plugins, which means running etcd plus APISIX yourself and writing the routing and plugin config by hand. My problem is "one person and a handful of agents on one machine", and deploying a gateway built for production traffic to solve it is a cannon for a mosquito.

## Still broken

- **No standard.** `responses` / `chat completions` / `messages` — the gateway just keeps translating forever. This isn't any one project's bug; it's an industry-wide tax on incompatibility (I ranted about it in {% post_link API-互不兼容-OpenAI和Anthropic的打印机耗材垄断策略 [the printer-ink monopoly post] %}).
- **Config is still scattered.** The gateway unified the endpoint, not the model list. The real fix is one source of truth generating every client's config — I haven't built it.
- **Zero cost observability.** The gateway routes; it doesn't tell me where the money goes.
- **No health checks on the chain.** With four layers, any one of them dying looks identical: "the agent stopped".

## Wrap-up

In the end I didn't actually achieve "not depending on an API" — I depend on *more* vendors than when I started. But that was never the real goal, and the real goal I did hit: the entry layer is mine now. Swapping vendors is one line in a config instead of rewriting five client setups.

The cost is three or four extra components and a debugging chain nobody enjoys reading about. Whether that trade was worth it, I still can't say — there's a decent chance next year's sequel is about me tearing it all down.
