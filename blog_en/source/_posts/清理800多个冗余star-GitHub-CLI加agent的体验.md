---
title: Cleaning Up 800+ Stale Stars — A GitHub CLI + Agent Experience
categories: Others
date: 2026-08-05 22:00:00
tags: ['GitHub', 'GitHub CLI', 'gh', 'AI Agent', 'star', 'Multica', 'cleanup']
---

The way I star GitHub repositories is pretty much the same as how I like and collect videos on bilibili: I just click them on a whim. Whenever I see a project that looks even remotely interesting, my brain always goes "might be useful someday" or "I'll take a look later" — and my finger clicks the button before my brain finishes the thought.

That's how, over the years, my star count quietly climbed past 800. To be honest, I have no idea what those 800+ projects actually are — most of them were starred because they seemed useful back then, and then never opened again.

<!-- more -->

## I'd thought about cleaning up before, but it was too much hassle

It's not that I never wanted to clean up. I've tried some star-management tools before — the kind that pull out your star list and group it by tag. But they were still a pain: after the sorting, I still had to look at every single one and click through it by hand. Going through hundreds of stars one manual operation at a time would wear my finger out. So the cleanup kept getting postponed, all the way until today.

## With an agent, everything simplifies

This time I simply delegated the job to an agent on [Multica](https://multica.ai/): pull out all my stars, match them against my existing 14 star lists, suggest tags for the untagged ones, and finally produce a complete table.

The agent's workflow was clear:

1. Pulled all 808 of my stars via the GitHub API (`gh api`)
2. Matched them against my 14 existing star lists (AI Or ML / Self-Host / Bioinformatics / AppDev / DataScience / We-media / etc.)
3. Kept the 309 already-tagged stars as-is, then judged the remaining 499 untagged projects one by one by their descriptions, suggesting tags that all fell into the existing 14 categories
4. Output everything as a markdown table for me to review

After I got the table, my part was simple: **keep what I want, delete what I don't**. Then I told the agent to batch-unstar everything not on the list.

Results:
- 335 stars were unstarred in total, bringing the count from 808 down to 473, exactly matching my keep-list
- For the kept projects, the 174 that needed tagging — 274 tag entries in total — were also written back to the corresponding star lists by the agent in one go

The operations I needed to do were very few: open coidum, review the table, delete the projects I didn't want. Everything else — fetching, matching, batch-unstarring, writing back tags — was done by the agent.

## Some thoughts

Those dedicated star-management tools are probably going to lose a chunk of their market. The core problem they solve — organizing scattered stars into categories and helping me bulk-clean them — a general-purpose agent can now do better. And the agent's advantage isn't just automation: it can intelligently judge which category a project belongs to based on its description, call the GitHub API directly for batch operations, and even caught a subtle detail for me along the way (that writing back to star lists requires the GraphQL API plus the `user` scope).

I'll keep clicking stars on a whim, of course — but at least the cleanup part won't give me a headache anymore.
