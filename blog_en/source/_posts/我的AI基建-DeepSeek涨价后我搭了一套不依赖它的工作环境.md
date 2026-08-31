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
  - Multica
  - MCP
  - Self-hosted
---

My day-to-day work is heavily multi-model and multi-agent now: several agents running in parallel across different projects, backed by a whole shelf of model vendors. Running a lot of work like this — reliably *and* cheaply — can't be held together with ad-hoc config edits. So after DeepSeek's price hike, instead of spending more effort on "how do I spend less", I flipped the question: **I need my own infrastructure, so that any single vendor failing (price hike, slowdown, protocol change, bad afternoon) costs me one leg, not the whole house.**

<!-- more -->

## What this infrastructure has to provide

Once thought through, the requirement splits into four pieces, each answering a different class of real incidents:

**1. Agent runtime.** My work laptop loses power and network — and when it does, half-finished tasks die with it. Agents have been my main workforce for a long stretch now; parking that workforce on a laptop that can trip a breaker is untenable. Tasks need to live on a board, get claimed by a runtime, and execute in a container independent of my machine — if the machine dies, the task stays in the queue. For this layer I use {% post_link Multica-AI原生的任务管理平台 [Multica] %} (my assessment of the idea is in {% post_link 开源平替Claude_Managed_Agents的Multica到底是什么 [this post] %}; the container setup is written up {% post_link 把Claude-Science装进Docker容器 [here] %}).

**2. Stable, acceptably-priced model vendors.** No single basket. The price hike is one thing; the everyday problem is "this vendor is lagging today". Every token-plan vendor has its bad days and rate-limit days, so: multiple vendors, multiple accounts, use whoever is healthy. The selection logic and the cost math (including each vendor's "sweetest tier" paradox) are in {% post_link Omniroute-DeepSeek涨价后-一个Key管起所有Token-Plan [the previous post] %}; the price-war backdrop in {% post_link DeepSeek-永久降价-大模型价格战的新玩法 [this one] %}; my usage scale {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [here] %}.

**3. The agent protocol layer.** The most underrated piece. Third-party token vendors sell you inference and **do not care about protocols**: Codex speaks OpenAI's `responses`, Claude Code speaks Anthropic's `messages`, opencode wants both. Dialects don't match, so the translation has to happen in a layer I own.

**4. An MCP gateway (not built yet, but planned).** The more tools an agent can reach, the better — but tools have management costs too: dozens of MCP servers' configs, auth, and versions scattered across agent configs will spiral eventually. I want another layer that aggregates and governs tools. There are two routes for tool reuse, skills and MCP; I lean toward MCP — it's protocol-level reuse, naturally portable across agents, while skills feel like one agent's private assets.

Put together, the trend is clear: my AI infrastructure is not "one gateway project" — **runtime, vendors, protocol, and tools each need their own owner**, each layer mine, with vendors as the outermost replaceable skin.

## The pitfalls in each layer

### Runtime: the environment is where the pain is

The pitfalls here are mostly not Multica's — they're "what's installed in the container and how". Agents run in containers, so the environment must be reproducible; I manage it with {% post_link 终于用pixi搞定了用于自动构建的recipe和action [pixi] %}, with the full build-recipe story in that post. The other recurring one is the permission model: agents with built-in sandboxing like Reasonix assume a full OS and need sandbox features disabled to run inside a container. Every new agent entering the container pays this tuition again.

### Vendors: old pitfalls unfixed, new ones keep coming

omniroute's three pitfalls (search interception not firing, brand-new vendors not syncing quota, heavy concurrent tasks disabled by default) are already covered in {% post_link Omniroute-DeepSeek涨价后-一个Key管起所有Token-Plan [the previous post] %}, so I won't repeat them. The real lesson of this layer: **multi-basket costs more to maintain than expected** — more baskets means more glue between baskets, and the glue itself becomes a new single point of failure.

### The protocol layer: the one a weekend couldn't fix

This is where I spent the most time and where this article has its most pitfalls. In detail:

**Generation 1: direct to the official API.** Container config pointed straight at DeepSeek's official API — Claude Code via Anthropic mode, Codex via OpenAI mode, plus a separately configured vision fallback model. Functionally fine. But base_url was scattered across every agent's config, and each new model in the gateway meant editing every client one by one. At one point I did a big flip: all traffic back to the official API, default model switched to a vision experimental build (keeping the stable-looking name in the UI). The day I counted how many places that touched — that's when I decided not to live like that anymore.

**Generation 2 (omniroute) solved "which vendor" but not "which protocol".** It unified vendors into one entry point and one key — genuinely nice — but left the protocol debt untouched: it assumes client and upstream speak the same dialect, and mine don't.

**Generation 3 died on `responses`.** I moved to new-api expecting more thorough protocol unification; it turns out upstream new-api doesn't support OpenAI's `responses` protocol — and Codex speaks precisely that. Configure the direct address locally, and the request simply isn't supported. No config workaround: the two protocol stacks don't overlap. (This incompatibility is an industry-wide tax — I ranted about it in {% post_link API-互不兼容-OpenAI和Anthropic的打印机耗材垄断策略 [the printer-ink monopoly post] %}.)

**Generation 4: add CLIProxyAPI purely to translate.**

```
Codex → CLIProxyAPI(:8317) → new-api → vendors
```

Its default port is 8317. On the Codex side you define a custom provider in `~/.codex/config.toml`, point `base_url` at it, and set `wire_api = "responses"` — that last one is Codex's protocol switch and the entire reason this layer exists. [Docs here](https://help.router-for.me/introduction/quick-start.html). My full chain now:

```
agents in the container (Codex / Claude Code / opencode / Reasonix / dsh / codebuddy)
   → CLIProxyAPI (protocol translation)
   → omniroute or new-api (routing and quota)
   → vendors
```

It works. The price: **every layer has its own config and its own defaults, so fault localization became "cut the chain layer by layer" elimination**. A good share of my debugging time wasn't fixing problems — it was establishing *which layer* was broken.

**bifrost and apisix: evaluated, never deployed.** Full disclosure — I never actually ran either, so there are no invented war stories here, only why I didn't ship them. bifrost is a high-performance Go LLM gateway, and performance is exactly what someone burning my volume should care about; but the capability I actually need — one entry feeding both OpenAI-style and Anthropic-style clients — is precisely where gateways differ most, and you can't tell from a README. You need your own traffic on it. When I evaluated it, omniroute was already carrying daily load, so the evaluation stayed an evaluation. apisix is a different species: a general-purpose API gateway with AI capabilities as plugins, meaning you run etcd plus APISIX yourself and hand-write routing and plugin config. My scenario is one person and a handful of agents on one machine; a gateway designed for production traffic is a cannon for a mosquito.

### MCP gateway: not built, so no pitfalls

I won't invent what doesn't exist. This layer sits last on purpose: **with the protocol layer still unstable, moving the tool layer first is wasted effort.** Once the gateway settles, I'll decide between building my own aggregator and adopting an off-the-shelf MCP gateway.

## Infrastructure is hard — and it's the only way to keep work moving

The protocol layer above: I spent one weekend on it and didn't solve it. Half of a second weekend more, and the chain barely worked — then began "which layer is broken" debugging. Even now it hasn't reached a finish line I'd call satisfying. This kind of investment is a luxury in daily work: it produces zero business code, pure paving for future work.

But run the other ledger. I computed agent ROI seriously in {% post_link 对Agent收益的思考 [this post] %}: an agent's value hinges on running *uninterrupted*. One vendor slowdown, one price change, one protocol mismatch — and the whole pipeline stops. **The value of infrastructure shows up precisely in the days you don't notice it**: nothing happens, tasks finish, and you forget you were debugging a gateway last week.

The current setup is still far from elegant: a four-layer chain, no single source of truth for config, zero cost observability, no health checks. But "swap vendors tomorrow" went from "rewrite every client's config" to "change one line, maybe patch one protocol translation".

Infrastructure that one weekend couldn't finish — next weekend, keep building.
