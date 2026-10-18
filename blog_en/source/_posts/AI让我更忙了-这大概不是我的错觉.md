---
title: AI Made Me Busier, and That's Probably Not Just Me
categories: AI
date: 2026-10-18 21:00:00
tags: ['AI', 'Agent', 'Anthropic', 'Workflow', 'Thoughts']
---

For the past couple of months I've been digging through my own numbers: commit counts, token consumption, how my tasks spread across the day. The conclusion is unambiguous — my workload really did go up, and not by a little. But I kept sitting on one doubt: is this a personal problem (bad task management, inability to say no, failing at managing my energy), or is it a broader shift that's happening to everyone?

On October 16, The Anthropic Institute published [How employment is changing at Anthropic](https://www.anthropic.com/research/employment-at-anthropic) (by Nathan Wilmers, Regina Catipon, Maxim Massenkoff and Peter McCrory). Using the company's own employment records, usage telemetry, output metrics, employee surveys, calendars and hiring scorecards, it describes how AI changed the shape of work inside one company. Reading it, most paragraphs had a matching line somewhere in my own work logs.

<!-- more -->

## Three Things I Felt First

### 1. The workload really did grow

I once [pulled two years of commits and analyzed their distribution](/en/2026/10/17/我给AI排了24小时班-把自己排成了值班表/). Same window (January 1 to August 31, 243 days):

| Item | 2025 | 2026 |
|---|---|---|
| Commits | 324 | 1013 |
| PRs | 150 | 345 |
| Days with activity | 111 | 171 |

August alone is even more striking: 23 commits and 10 PRs for the whole month last year, versus 167 and 84 this year. And the second jump in the curve lands in May 2026 — exactly when AI started doing a large share of the work.

Token consumption tells the same story: I recorded [38M to 500M in three months](/en/2026/08/05/三个月Token消耗暴涨13倍-我的AI工作流演变记录/), 883M in total.

In one sentence: **AI didn't make me do less work. It made me capable of doing more, and then I actually did more.**

### 2. My brain isn't enough

In [Multiple Devices, Multiple Platforms, Multiple Projects — Too Bad I Don't Have Multiple Brains](/en/2026/06/14/多设备-多平台-多项目-可惜我没有多个脑子/) I summed up my own state: **async works, multi-threading doesn't.** Agents can efficiently preprocess a lot of the groundwork, but every task still ends up queued behind me — I can't review and handle several different tasks at the same time.

### 3. I became the bottleneck

In that same stats post I wrote: my review-and-judgment speed can't keep up with their output. That's why the 24-hour distribution didn't get noticeably longer in 2026 — the actual writing and research isn't mine anymore — but the early-morning bars did get taller.

## Then I Hit This Passage in the Report

The report has a counterintuitive finding: **in interactive sessions, roughly 80% of the time the agent is waiting on human input**, and frontline engineers describe people as "bottlenecking the agents".

I paused there. When I wrote "my review speed can't keep up with their output", I assumed it was just me — insufficient ability, badly decomposed tasks. And then a frontier lab went and measured the same thing with telemetry, and landed on almost the same sentence as my complaint.

Following that thread, a few more numbers:

- **Meta-work is ballooning**: sessions spent writing specs for agents, having one agent review another, or building monitoring dashboards went from **10%** in autumn 2025 to **22%** in spring 2026.
- **Standing instructions pile up**: the number of standing instructions written for agents went from a few hundred in July 2025 to **8,500** a year later.
- **Meetings went up, not down**: median employee meeting time is **+20%** year over year, and **+16%** for technical ICs, which the report explains as trust and accountability still needing face-to-face coordination.

These hit close to home, because that's basically what my 2026 looked like: my time moved from "writing code" to "writing rules, writing specs, writing acceptance criteria, and cleaning up after agents".

| What the report observed | What I lived |
|---|---|
| ~80% of interactive session time is agents waiting on people | "My review speed can't keep up with their output" |
| Meta-work sessions 10% → 22% | Time shifted from coding to spec-writing, dispatching and review |
| Standing instructions: hundreds → 8,500 | Workspace rules and prompt templates keep accumulating |
| Commits per engineer rose to 5× in six months | My commits went 324 → 1013 from Jan–Aug |
| Median meeting time +20% | Explaining "what the AI actually finished" costs more words |

## "6×" Is Not 6× Output

The headline-friendly number: by August 2026, the interactive agents commanded by an individual employee logged **six times** that employee's own hours. But read the definition:

- agent hours are counted as **runtime** (including reasoning and tool calls, excluding any single wait longer than 30 minutes);
- human hours are estimated as **headcount × 40 hours per week**;
- total agent hours first exceeded total human hours in March 2026.

So it measures **machine-running time**, not 6× output, and certainly not 6× headcount replacement. Don't let the headline run away with you.

The genuinely interesting item sits in the efficiency section: comparing human-written and agent-written code for the **same employee, same month, same directory**, the share **deleted or reverted within 90 days showed no significant difference**. So "everything AI writes is garbage" doesn't hold either.

Put the two together and it gets a little paradoxical: AI is fast and not obviously worse, so why am I more tired?

The report's answer is pretty plain: **the bottleneck isn't generation, it's judgment. And judgment can't be parallelized, outsourced, or cached.**

## The Part About Junior Roles

The sharpest structural finding:

- Share of **junior employees** in technical orgs: 21% in early 2025 → 14%;
- Share of **junior technical hiring**: 27% → 9%, with postings aimed at juniors dropping especially steeply after December 2025;
- Technical orgs went from 61% to 56% of headcount, but technical hiring **in absolute terms is still growing** — sales, finance and recruiting just grew faster.

"Total volume is still growing" is no comfort to the segment being squeezed out structurally. For those already aboard, the report supplies the other half: 27% of employees say Claude clearly made them more productive, 20% say their work went "from doing it myself to managing Claude", and employees with a year of tenure are 50% more likely per month to move between departments in 2026 than in 2025.

## But Don't Take It as Proof

The authors list the limitations themselves, and I think they matter more than the numbers:

- **Not causal**: tool adoption isn't random — people about to start a big project are more likely to arm themselves first, so there's selection bias.
- **One company, extremely unusual**: Anthropic tripled in size within a year and its revenue per employee is 2–4× the peak of Meta's and Google's high-growth eras. Its practices can't be extrapolated to ordinary companies.
- **It doesn't address people who were never hired**: internal redeployment policies do nothing for those who never got an offer.

So I won't write "this will happen to every worker". The careful version: this is one company opening up its internal data, and the mechanism it reveals — **people as the bottleneck, meta-work ballooning, coordination costs rising** — matches my own records. How strongly that mechanism plays out elsewhere depends on the kind of work and the depth of AI adoption at that company.

But it does answer part of my opening question: if even Anthropic's employees burn half their time queueing agents and writing agent instructions, then my exhaustion probably isn't simply me failing to manage myself.

## So What Can We Do

No grand narrative — just the parts a worker can still move:

**1. Admit people are the bottleneck, then staff around the bottleneck.**
If the constraint is my judgment, opening three or four agents in parallel won't raise output — it only lengthens the queue. That 80% figure is essentially saying agent throughput is already capped by human review speed.

**2. Treat meta-work as real work.**
Writing standing instructions for your agents, defining acceptance criteria, keeping checklists — these aren't "preparation", they're a main part of the job now. Those 8,500 standing instructions are exactly this point: that stuff is an asset, worth accumulating.

**3. Protect your judgment and your context.**
In hiring scorecards, technical roles increasingly weigh **autonomy and judgment**, while domain expertise and execution weigh less. That signal cuts both ways for an individual: pure execution seats are shrinking, but if your value is "knowing what should be done and what shouldn't be accepted", you're worth more, not less.

**4. Don't mistake tokens and commits for achievement.**
My 13× tokens and 3× commits mean "I used AI to do more things", not "I produced 3× the value". People who conflate the two end up fooling both their boss and themselves.

**5. Watch out for "so-so automation".**
The report warns about automation that replaces people without delivering much real productivity. At the individual level: if an AI finishes something but I still spend the same time verifying it, the work didn't disappear — it just changed from "do it" to "check it".

## At the End

The report's three employer actions are broad-based equity, job early-warning systems, and prioritizing internal redeployment. All of that lives at the company level, beyond anything I can influence as an employee.

What I can influence is probably just that line from an older post: **instead of complaining, find a way to adapt.**

As for how well that adaptation goes — let me at least finish this blog post first.
