---
title: The Benefits and Imperfections of Claude Science
date: 2026-07-15 22:00:00
tags:
  - Claude Science
  - AI
  - Review
  - Bioinformatics
categories: AI
---

After previously tinkering with putting Claude Science into a Docker container and using it for a while, I want to talk about its benefits and the areas where it still falls short.

<!-- more -->

## Benefits: Rich Ecosystem of Biology Plugins

What I'm most satisfied with about Claude Science is how comprehensive its biology-related plugins are. No special configuration is needed—you can start searching for data, setting up download environments, and performing analyses right away.

It comes with built-in scientific MCP servers such as biorxiv and chembl, allowing you to pull data directly from public databases. For people working in bioinformatics, this means saving a huge amount of time on environment setup. In the past, just installing software, configuring the environment, and finding data for an analysis could take half a day. Now, with Claude Science, you can get a workflow running with just a few sentences.

Its support for commonly used bioinformatics tools is also decent. I tried using it for simple single-cell analyses, and with the built-in plugins, I basically didn't need to prepare any code myself—just describing the requirements clearly was enough. For quickly validating ideas, the experience is indeed very good.

## Imperfections: Reproducibility Falls Short of Its Name

Claude Science adopts a cell-based editing and execution approach similar to Jupyter notebooks. The advantage is that users can view variables in the runtime environment in real time, providing a strong sense of interactivity—this is definitely commendable.

However, the problem is that **after multiple rounds of conversation, the code retained in the interface is not the code that actually produced the results**.

Anyone who has used Jupyter knows the classic scenario: you repeatedly modify and run the same cell, and eventually the cell displays the last edited code, but the results from earlier runs may correspond to completely different code. Claude Science has a similar issue—the conversation history retains the final state of the code, not the version of the code that "actually produced the results."

This leads to a problem: **the reproducibility it claims at launch is compromised in actual use**. When you look back at your previous analyses, it's hard to confirm whether the code you see matches the results at that time. If you accidentally close the session and reopen it, reproducing the previous results becomes even more difficult.

## Summary

Claude Science does a good job of lowering the barrier to bioinformatics analysis. Its out-of-the-box plugins make some routine analyses very efficient. However, its shortcomings in code version management make its analyses and results unreliable.

At its core, it's better suited as an exploratory analysis tool—quickly validating ideas and doing preliminary data processing, where the experience is great. But for research scenarios that require rigorous recording of the analysis process and truly reproducible results, there is still considerable room for improvement. A good tool should not only be convenient to use, but also trustworthy. I hope Anthropic can address this in the future.
