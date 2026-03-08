---
title: Configuring Multi-Agent Workflows with OpenCode
categories: Coding
date: 2026-03-08 16:17:24
keywords: [opencode, subagent, skill, multi-agent workflow, AI workflow, automation]
---

Last year, I explored AI topics purely out of personal interest. This time, it's under pressure from my boss... Even though I'm full of resentment, I still need to write something down to accumulate experience...

<!-- more -->

The task I need to accomplish is to create a workflow using AI tools. This workflow should enable AI to automatically perform web searches based on specific keywords, summarize the retrieved content, and recommend topics for our company's official WeChat account. This is a relatively simple workflow, but it involves two critical issues:

1. **Automatic verification of information authenticity**: The content will be published on the company's official account, so it must be verified and fact-checked. However, no matter which model is used, hallucinations cannot be 100% avoided. If we cannot impose maximum restrictions, the time cost of checking whether the AI is hallucinating would be higher than manually selecting topics, which would defeat the purpose of this workflow.

2. **Performance limitations of workflow tool agents**: Although the model handles the key processing, after using several coding agents, I've realized that the design of the agent program itself significantly affects the effectiveness and efficiency of the final output. Tools like n8n and dify allow setting up agent nodes, but these agents are simple ones developed or configured by developers. Their usability is certainly far lower than well-designed systems like OpenCode, Claude Code, or OpenClaw, and their effectiveness is likely greatly reduced.

Therefore, to quickly produce a demo within one day (yes, a requirement casually mentioned in a weekend group chat, with only one actual working day available... in the past, even understanding the basic technology would take longer), I decided to try the Subagent mechanism built into OpenCode.

## Detailed Explanation of the Subagent Mechanism

### What is a Subagent?

OpenCode has two types of agents: **Primary** and **Subagent**. The primary agent is the one users directly interact with, while subagents are specialized assistants that the agent can call to perform specific tasks.

One important purpose of subagents is context isolation. Even though today's model context windows have expanded several times compared to the past, a single conversation's context is still insufficient to hold all the code and execution details of a task. Therefore, splitting related but not directly relevant content into subagent conversations, and keeping only the work summaries and results of subagents in the main agent's conversation, helps the agent stay on track and avoid deviation.

### Built-in Subagents

OpenCode comes with two built-in subagents:

1. **General**: A general-purpose agent for researching complex problems and executing multi-step tasks. It has full access to tools (except todo), can modify files when needed, and can be used to run multiple work units in parallel.

2. **Explore**: A fast read-only agent for exploring the codebase. It cannot modify files. Use it when you need to quickly find files by pattern, search for keywords in code, or answer questions about the codebase.

### Primary Agents

In addition to subagents, OpenCode also includes two primary agents. In the TUI interface, you can switch between them by pressing Tab:

- **Build**: The default primary agent, with all tools enabled. Used for development work requiring full access to file operations and system commands.
- **Plan**: A restricted agent that, by default, prohibits file editing and bash commands. Used for analysis and planning without making any changes.

### Custom Agents

Both types of agents (primary and subagent) support customization. You can create custom agents in two ways:

#### Method 1: Using Markdown Files

Place custom agents in the specified directories:
- Global: `~/.config/opencode/agents/`
- Project-level: `.opencode/agents/`

The filename becomes the agent name. For example, creating `review.md` creates an agent named `review`:

```markdown
---
description: Reviews code for best practices and potential issues
mode: subagent
model: anthropic/claude-sonnet-4-20250514
temperature: 0.1
tools:
  write: false
  edit: false
  bash: false
---

You are in code review mode. Focus on:
- Code quality and best practices
- Potential bugs and edge cases
- Performance implications
- Security considerations
```

#### Method 2: Using opencode.json Configuration

```json
{
  "agent": {
    "review": {
      "description": "Reviews code for best practices",
      "mode": "subagent",
      "model": "anthropic/claude-sonnet-4-20250514",
      "prompt": "{file:./prompts/code-review.txt}",
      "tools": {
        "write": false,
        "edit": false
      }
    }
  }
}
```

#### Key Configuration Items for Custom Agents

- **description** (required): A functional description of the agent for the AI to decide when to call it.
- **mode**: Agent mode, options: `primary`, `subagent`, or `all`.
- **model**: Specifies the model to use.
- **temperature**: Controls output randomness (0-1).
- **tools**: Controls the tools available to the agent.
- **prompt**: Path to a custom system prompt file.
- **hidden**: Set to `true` to hide the subagent from the @ autocomplete menu.

## Detailed Explanation of the Skill Mechanism

### What is a Skill?

A Skill is OpenCode's reusable instruction mechanism, allowing you to define reusable behaviors through `SKILL.md` files. Skills are loaded on demand via the native `skill` tool—agents can view available skills and load their full content when needed.

### Skill Placement

OpenCode searches the following locations to discover skills:

- Project configuration: `.opencode/skills/<name>/SKILL.md`
- Global configuration: `~/.config/opencode/skills/<name>/SKILL.md`
- Project Claude-compatible: `.claude/skills/<name>/SKILL.md`
- Global Claude-compatible: `~/.claude/skills/<name>/SKILL.md`

### Skill Structure

Each `SKILL.md` must start with YAML frontmatter containing the following required fields:

- `name`: Skill name (1-64 characters, lowercase letters and numbers only, may be separated by a single hyphen).
- `description`: Skill description (1-1024 characters).

Example:

```markdown
---
name: web-search
description: Skill for web searching
license: MIT
---

## What I do
- Perform web searches based on user needs
- Retrieve and summarize search results

## When to use me
Use this skill when users need to query real-time information or the latest news.
```

### How to Use Skills

During configuration, you can decide which skills the agent needs to preload. Loading essentially means placing the content of SKILL.md into the context.

```
skill({ name: "web-search" })
```

Additionally, most skills describe when they should be invoked, so if the conversation prompt matches the conditions for tool invocation, the skill will be automatically loaded (unless explicitly denied in the definition).

## Practical Application

In this project, the steps can be divided into three parts:

- Information retrieval and collection
- Topic selection
- Information verification

Then, a primary agent is used to control the aforementioned workflow.

## Problems Encountered

### Source Link Distortion

This problem was more severe than I anticipated.

I had the agent first retrieve information via Agent-Browser, then generate content based on the retrieved information. These parts actually worked quite well, and the output wasn't too distorted.

However, the seemingly simple task of saving source information links was something that Minimax-M2.5, used during initial development, simply couldn't do correctly. Nine out of ten saved links were wrong. Even with the temperature set to the minimum and additional prompts added, it couldn't preserve the real links.

When cross-checking information, the situation was even more absurd. When saving links, it would most likely modify the singular/plural forms of words in the links and some code strings. For example:

- video > videos
- a-s001 > b-a273

After struggling for a long time, switching to GLM resolved this issue...

### Unstable Agent Actions

This problem affects various aspects of the generation process.

For the primary control agent, when it schedules multiple subagents, its control stability isn't very flexible. For instance, my workflow has three steps. Sometimes I want to execute only part of them, but if my input instruction is just a small question within the conversation prompt, the model might still run the entire process based on the longer predefined prompt.

For subagents, once they leave the previous context or if the implementation method isn't explicitly specified, each implementation could be completely different. For example, writing JavaScript this time, Python next time, and directly executing Shell commands another time. The whole process resembles the early days of text-to-image models with repeated "drawing lots"...

## Lessons Learned

- Different models have different strengths and should be used in combination.
    + GLM: Accurate instruction execution, but output is very slow.
    + Minimax: Unstable for non-coding instructions, but good at writing code, with fast output.
    + Kimi: Moderate in both instruction execution and coding, but prone to inexplicable infinite loops, making it unsuitable for control agents.

- For steps with high determinism, using agents is less efficient than directly generating code and skills, then having the agent call them.

- Workflows still have their value. Compared to multi-agent systems based on strong agents, workflows offer better stability at the cost of some flexibility. The choice should be made based on the specific task.
