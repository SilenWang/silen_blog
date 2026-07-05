---
title: Putting Claude Science into a Docker Container
date: 2026-07-05 23:00:00
tags:
  - Docker
  - Claude Science
  - API
  - AI
categories: AI
---

I recently worked on a small project: putting the Linux version of Claude Science into a Docker container, paired with an API Bridge so it can run through third-party API backends without needing to hold an Anthropic API Key directly.

The process was more eventful than expected—also quite interesting. Here's a record of it.

<!-- more -->

## Project Architecture

The overall architecture is actually not that complicated:

```
┌─────────────────────────────────────────┐
│              Docker Container            │
│                                          │
│  ┌──────────────┐    ┌────────────────┐  │
│  │  API Bridge  │◄───│ Claude Science  │  │
│  │  (:9876)     │    │  (:9981)        │  │
│  └──────┬───────┘    └────────────────┘  │
│         │                                │
└─────────┼────────────────────────────────┘
          │
     ┌────┴────────────┐
     │ Third-Party API │
     │ (SiliconFlow,   │
     │  Moonshot, etc) │
     └─────────────────┘
```

Claude Science runs via the operon daemon, which sends Anthropic-format requests to the [API Bridge](https://github.com/SilenWang/Claude-Science-Container). The API Bridge converts requests to OpenAI-compatible format and forwards them to the backend third-party API.

## Pitfalls During Development

### 1. Access Address Restriction

Claude Science provides a `--host` parameter to specify accessible addresses, but in practice, you cannot access the software interface from any IP other than localhost. So either you use port forwarding for local use, or you need to put the tunneling service and this container on the same machine—otherwise it simply won't work...

### 2. HTTPS URL Restriction

Originally I followed the principle of one service per container, separating the API Bridge and Claude Science. But I discovered that Claude Science restricts request URLs to HTTPS or localhost, so I had to fall back to a single container...

### 3. Protocol Conversion Limitations

When the API Bridge does the Anthropic -> OpenAI format conversion, some Anthropic-specific features are lost. The most obvious are **tool calls** and **web_search**. OpenAI-compatible endpoints don't have corresponding implementations, so these features become unavailable.

That said, DeepSeek provides an Anthropic-compatible direct API endpoint. If you go through this channel, you can bypass the conversion and keep all features intact. But the trade-off is a narrower model selection—you're limited to DeepSeek's official model API.

### 4. MCP Directory Connector Issues

When entering Claude Science's settings interface, it shows "Directory connectors unavailable". At first I thought it was a configuration issue, but after checking the logs, I found `claudeAiFetch: 401 and refresh failed`—it was trying to sync the MCP server list from claude.ai, but since I wasn't logged into an Anthropic account, it failed.

However, this doesn't actually affect core functionality. Claude Science's built-in scientific MCP servers (biorxiv, chembl, etc.) still connect and work fine. It's one of those "error reported but doesn't affect usage" situations.

This reminded me of using snakemake back in the day—some rules' shadow functionality would show warnings in certain versions without affecting execution. I found those semi-broken features annoying at the time, but I also learned how to quickly determine which warnings to ignore. This time with the MCP directory error, one glance at the logs told me it was an auth issue, not a functional problem, so I didn't dwell on it.

### 5. Unexpected Task Termination

When running tasks with DeepSeek Flash, tasks would sometimes terminate abruptly mid-execution, but DeepSeek Pro was much better. This issue was hard to pin down—it's neither a configuration problem nor a network issue; it seems more like behavioral differences in the model itself. I can only assume the replacement model isn't quite smart enough yet...

## Project Repository

The code is on GitHub: [Claude-Science-Container](https://github.com/SilenWang/Claude-Science-Container)

Feel free to give it a try if you're interested. Open an issue if you run into problems.

## A Small Rant

I wrote a blog post before complaining about API incompatibility—how OpenAI and Anthropic APIs aren't compatible, and switching models requires rewriting a bunch of code. I compared it to printer manufacturers' consumables monopoly tactics.

Little did I expect that complaint would quickly turn into a real problem I had to deal with. Claude Science, being an Anthropic product, naturally uses Anthropic-style API calls. But I wanted to use third-party API backends (like SiliconFlow, Moonshot, etc.), which only support OpenAI-compatible formats. This means a translation layer is needed in the middle to convert Anthropic-style requests to OpenAI-compatible format. Fortunately, the open-source community came through—a solution appeared overnight, and all I had to do was be a package integrator to make it work.
