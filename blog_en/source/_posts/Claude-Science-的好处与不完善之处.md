---
title: The Benefits and Limitations of Claude Science
date: 2026-07-15 22:00:00
tags:
  - Claude Science
  - AI
  - Bioinformatics
categories: AI
---

After previously tinkering with putting Claude Science into a Docker container, I've been using it relatively intensively over the past two weeks (connected to DeepSeek, not the original Anthropic version, of course), so I want to talk about its advantages and limitations.

<!-- more -->

## Advantages: A Well-Equipped Biomedical Workbench

What I'm most satisfied with about Claude Science is how comprehensive its biology-related plugins are. Without any special configuration, you can start searching for data, setting up download environments, and performing analyses.

It comes with built-in database connectors (MCP servers) for the most commonly used resources in biomedical research. Without any special setup, you can search for literature on PubMed, find validation datasets on GEO, and use Conda to deploy the tools needed for analysis. For anyone looking to validate ideas using existing data and classic methods, Claude Science can save a significant amount of time. In the past, just installing software, configuring the environment, and finding data for an analysis could take half a day. Now, with Claude Science, you can get a basic workflow running with just a few sentences.

I tried using it for simple single-cell analyses, and with the built-in plugins, I basically didn't need to prepare any code myself—just describing the requirements clearly was enough. For quickly validating ideas, the experience is indeed very good.

Additionally, the built-in Reviewer mechanism is very useful. When Review is enabled, the output text, the internal logic of the results, etc., are all automatically checked. If an issue is found, it will automatically redo the work—at the cost of doubling the token usage, of course.

## Limitations: Not a Reliable Analysis Expert

Just like when I first started Vibe Coding, after the initial burst of joy, you gradually start to notice various shortcomings...

I should note upfront that I'm using a version accessed via a community project connected to DeepSeek, so most of the issues in the limitations section are clearly compatibility problems. As for the remaining issues, it's hard to rule out that they are also caused by compatibility. After all, Anthropic's products are designed to work as a software-model package deal; once you decouple them, problems are inevitable.

### Incomplete Reproducibility

Claude Science adopts a cell-based editing and execution approach similar to Jupyter. The advantage is that users can view variables in the runtime environment in real time and verify the content being computed.

However, the problem is that, just like how people often execute Jupyter cells out of order, **after multiple rounds of conversation with Claude Science, the code retained in the notebook may not match the code that actually produced the results**.

This leads to a problem: **the reproducibility it claims at launch is compromised in actual use**. When you look back at your previous analyses, it's hard to confirm whether the code you see matches the results at that time. If you accidentally close the session and reopen it, reproducing the previous results becomes even more difficult.

So as a user, at certain stages it's still necessary to have Claude Science write the scripts, then run them manually to verify the results before interpreting them.

### GPU and Interactive Kernel Compatibility Issues

This is clearly a problem specific to container usage. Using NVIDIA's CUDA container, you can successfully pass through the GPU, but the `/proc` directory mount inside the container causes issues, making the interactive kernel unusable and breaking all components. I haven't found a workaround yet, so drug target prediction, structural prediction, and similar features are all currently unavailable.

### Limited Custom Skills

Adding custom skills and connectors requires online verification from Anthropic's servers. Naturally, this doesn't work with my third-party integration approach. The community project I'm using also doesn't include skill management features. Fortunately, the built-in functionality is rich enough that this isn't a major issue for now.

### Limited Image Understanding

Since DeepSeek is not yet a multimodal model, Claude Science's built-in image analysis capabilities don't work. I've already felt this limitation while analyzing spatial transcriptomics data.

### Plan Generation Failures

At the start of a research session, Claude Science, like a Coding Agent, parses the input and generates a plan. When using third-party models, there are clearly minor issues with format compliance, resulting in a certain probability that the generated plan doesn't meet the required format, causing the task to disconnect midway. However, this isn't a big problem—you can just continue the conversation directly.
