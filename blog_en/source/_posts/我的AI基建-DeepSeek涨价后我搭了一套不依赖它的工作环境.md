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
  - Claude-Science
  - Self-hosted
---

My day-to-day work is heavily multi-model and multi-agent now: several agents running in parallel across different projects, backed by a whole shelf of model vendors. Running a lot of work like this — reliably *and* cheaply — can't be held together with ad-hoc config edits. So after DeepSeek's price hike, instead of spending more effort on "how do I spend less", I flipped the question: **I need my own infrastructure, so that any single vendor failing (price hike, slowdown, protocol change, bad afternoon) costs me one leg, not the whole house.**

<!-- more -->

## What this infrastructure has to provide

Once thought through, the requirement splits into four pieces, each answering a different class of real incidents:

**1. Agent runtime.** My work laptop loses power and network — and when it does, half-finished tasks die with it. Agents have been my main workforce for a long stretch now; parking that workforce on a laptop that can trip a breaker is untenable. Tasks need to live on a board, get claimed by a runtime, and execute in a container independent of my machine — if the machine dies, the task stays in the queue. For this layer I use {% post_link Multica-AI原生的任务管理平台 [Multica] %} (my assessment of the idea is in {% post_link 开源平替Claude_Managed_Agents的Multica到底是什么 [this post] %}, and the daily-use gripes are in {% post_link AI工具使用体验与吐槽-Multica-MiniMax-Google-One-Gemini [this one] %}); the runtime itself is a container I built.

**2. Stable, acceptably-priced model vendors.** No single basket. The price hike is one thing; the everyday problem is "this vendor is lagging today". Every token-plan vendor has its bad days and rate-limit days, so: multiple vendors, multiple accounts, use whoever is healthy. The selection logic and the cost math (including each vendor's "sweetest tier" paradox) are in {% post_link Omniroute-DeepSeek涨价后-一个Key管起所有Token-Plan [the previous post] %}; the price-war backdrop in {% post_link DeepSeek-永久降价-大模型价格战的新玩法 [this one] %}; my usage scale {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [here] %}.

**3. The agent protocol layer.** The most underrated piece. Third-party token vendors sell you inference and **do not care about protocols**: Codex speaks OpenAI's `responses`, Claude Code speaks Anthropic's `messages`, opencode wants both. Dialects don't match, so the translation has to happen in a layer I own.

**4. An MCP gateway (not built yet, but planned).** The more tools an agent can reach, the better — but tools have management costs too: dozens of MCP servers' configs, auth, and versions scattered across agent configs will spiral eventually. I want another layer that aggregates and governs tools. There are two routes for tool reuse, skills and MCP; I lean toward MCP — it's protocol-level reuse, naturally portable across agents, while skills feel like one agent's private assets.

Put together, the trend is clear: my AI infrastructure is not "one gateway project" — **runtime, vendors, protocol, and tools each need their own owner**, each layer mine, with vendors as the outermost replaceable skin.

## One exception: the research workbench doesn't follow this rule

Those four points are the rules I set for myself — but there is one workbench I deliberately kept outside them: **Claude Science**.

It's {% post_link Claude-Science-的好处与不完善之处 [the workbench I use for research] %}; when I {% post_link 把Claude-Science装进Docker容器 [packed it into a Docker container] %}, a thin API Bridge swapped its backend to DeepSeek's official API in Anthropic-compatible mode, so I never had to hold an Anthropic key. It is **a separate thing** from the infrastructure above: not part of Multica, no shared runtime.

Why the exception? Because **research needs immediate interaction far more than other tasks**. Writing code, running batch jobs, tidying data — those can be "thrown onto the board, come back for the result"; an agent working for twenty minutes without me is fine. But looking at a plot, tweaking one parameter, seeing the result and deciding the next step — in that loop I must be present, and every added wait breaks the thread of thought. So I still use Claude Science directly instead of turning it into another agent task.

For the same reason, this workbench is the one corner of the whole setup that **connects to no third-party token vendor at all**. Not an oversight — a refusal. What I need from it is exactly the three things the official API happens to cover:

- **Multi-protocol**: Claude Science speaks the Anthropic dialect, and the official API happens to have a compatibility mode for it — the Bridge inside the container only adapts protocols, it touches neither routing nor quota.
- **Vision**: reading images is the daily bread of research work, plus I've set a dedicated vision-fallback environment variable, so when the main model can't see there's something underneath to catch it.
- **Web access**: literature lookups and page fetches depend on it, and third-party vendors mostly don't give you network tools at all.

This workbench is thus a living fossil of "generation 1": multi-protocol, vision, web — all **configured straight to the official API via environment variables in the container config**, with no routing, quota, or concurrency concerns, because the whole point is one layer of uncertainty less. The cost is real too: its config shares nothing with the main chain, so the same model name has to be maintained in two places.

Around it sits its own set of pitfalls, unrelated to the main line, just a few memorable ones. To run it in a container I wrote a Go launcher: multiple forwarding targets in a config file with interactive selection at startup, a configurable command, and it can run without the container too. Turning on GPU passthrough turned into its own disaster — pulling NVIDIA's official image killed every language kernel in the container, and the fix ran all the way from suspecting `privileged` being neutered through proc mounts to device mapping ({% post_link 容器内用opencode修复GPU直通 [this post] %}); GPU support ultimately goes through CDI mounts rather than baking CUDA into the image. Chinese text in Python/R plots inside the container rendered as tofu boxes and needed its own font fix. The wildest one: the forwarding program **got flagged as a virus and deleted by Windows** — a small Go binary I compiled myself, unsigned, with every antivirus its own imagination. It ended in researching code signing.

## The pitfalls in each layer

### Runtime: every new agent means paying the tuition again

Self-building the runtime container, the first day settled the obvious: agents don't run as root, so a separate user was created to execute them; the runtime's state directory and the agents' config directories are bind-mounted outside the container for persistence — otherwise every container rebuild means refilling the model-provider settings; and all dependencies go through {% post_link 终于用pixi搞定了用于自动构建的recipe和action [pixi] %} so the container stays reproducible.

Everything after that was tuition:

- **Installing one agent is a campaign of its own.** Today it's the Codex CLI and CodeBuddy; tomorrow it's Claude Code plus a CLI/TUI utility for switching keys, MCP servers, and providers between several assistants (the very existence of such a tool says something: client configs got messy enough to deserve a dedicated manager).
- **Sandbox features don't work in the container.** Agents that ship with built-in sandboxing assume a full OS is theirs; the first thing to do inside a container is often turning those features off, or they won't start.
- **Containers inside containers.** When an agent needs Docker, host-side "Docker-in-Docker" support has to be explicitly configured — it's not on by default.
- **Config has to actually externalise.** Originally the key and URL were hardcoded — one vendor only. Only later did it become a config file that can hold multiple keys and specify which vendor a given agent should use. Small step, but it's the client-side foundation of the whole multi-basket strategy.
- **Rules must apply to all agents at once.** I added a "write the least code" rule plugin for every agent in the container, configured in the image/init flow so each new container grows up with it. Don't underestimate this: agent behaviour rules enforced by each repo adding its own files will drift sooner or later.
- **Protocol pass-through bugs surface on this layer.** Running dsh with v4 pro: ordinary chat works, but entering thinking mode gets rejected upstream — `The reasoning_content in the thinking mode must be passed back to the API`. In thinking mode the reasoning content must be passed back to the API, and this request chain doesn't pass it back. *Same model, same vendor, different runtime — broken.* Gateways can't help here; you inspect request bodies layer by layer.

The real lesson of this layer: **a runtime isn't a program, it's a pile of conventions.** Every new agent reconnects the questions "where does its config live, what is it allowed to change, what permissions does it assume".

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
