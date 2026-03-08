---
title: 基于opencode的多agent工作流配置
categories: Coding
date: 2026-03-08 16:17:24
keywords: [opencode, subagent, skill, 多agent工作流, AI工作流, 自动化]
---

去年研究AI相关的内容都是出于爱好，这次是老板的压迫... 即使怨念深重，还是要来写点东西积累下经验...

<!-- more -->

我要完成的工作是，使用AI工具创见一个工作流，这个工作流需要让AI自动的根据特定关键词进行网络检索，然后小结检索到的内容，进行公众号选题推荐。这是个比较简单的工作流，不过其中涉及两个比较关键的问题：

1. 信息真实性自动核实：这些内容是要发到公司公众号上的，自然得是经过证实和核查的信息。然而不管是什么模型，幻觉是不可能百分百避免的，如果不能进行最大程度的限制，那核查AI有没有出幻觉的时间成本比人来选题还要高，就失去这个工作流的意义了。

2. 工作流工具的Agent性能限制：虽然完成关键处理的是模型，但是在使用过多个Coding Agent之后，我能感受到，Agent程序本身的设计还是非常影响最终产出结果的效果和效率的。我知道的n8n、dify这些，虽然可以设置Agent节点，但是这里的Agent，是纯靠开发人员自己开发或设置的简易Agent，它的可用性肯定远低于经过设计的opencode、calude code或openclaw，效果很可能大打折扣。

因此，为了快速在一天内能出个demo（是的，周末群里随口一说的需求，实际能给的有效工期就只有一天...搁过去了解基本技术都不够）,我决定用opencode内置的Subagent机制来试试。

## Subagent 机制详解

### 什么是 Subagent

OpenCode 中有两种类型的代理：**主代理（Primary）**和**子代理（Subagent）**。主代理是用户直接交互的agent，而subagent是agent可以调用来执行特定任务的专业助手。

子代理的重要初衷之一是隔离上下文，即使在模型上下问张读已经比过去扩展了数倍的今天，一次对话的上下文也不足以将一项工作的所有代码和执行内容存下来，因此，将有关联但是具体内容并不直接相关的内容，拆分到subagent的对话中，在agent的对话中只保留subagent的工作摘要和工作结果，有利于让agent在执行过程中保持目标，不要跑偏。

### 内置子代理

OpenCode 内置了两个子代理：

1. **General**：用于研究复杂问题和执行多步骤任务的通用代理，拥有完整的工具访问权限（todo 除外），可以在需要时修改文件，可用于并行运行多个工作单元。

2. **Explore**：用于探索代码库的快速只读代理，无法修改文件。当需要按模式快速查找文件、搜索代码中的关键字或回答有关代码库的问题时使用。

### 主代理

除了子代理，OpenCode 还内置了两个主代理，在Tui界面上，按Tab可以切换他们：

- **Build**：默认主代理，启用所有工具，用于需要完全访问文件操作和系统命令的开发工作
- **Plan**：受限代理，默认禁止文件编辑和 bash 命令，用于分析和规划而不进行任何更改

### 自定义代理

两种代理（主代理和子代理）都支持自定义。你可以通过以下两种方式创建自定义代理：

#### 方式一：使用 Markdown 文件

将自定义代理放在指定目录：
- 全局：`~/.config/opencode/agents/`
- 项目级：`.opencode/agents/`

文件名即为代理名称。例如，创建 `review.md` 会创建一个名为 `review` 的代理：

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

#### 方式二：使用 opencode.json 配置

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

#### 自定义代理的关键配置项

- **description**（必填）：代理的功能描述，供 AI 选择何时调用该代理
- **mode**：代理模式，可选 `primary`、`subagent` 或 `all`
- **model**：指定使用的模型
- **temperature**：控制输出随机性（0-1）
- **tools**：控制代理可用的工具
- **prompt**：自定义系统提示词文件路径
- **hidden**：设为 `true` 可将子代理从 @ 自动补全菜单中隐藏

## Skill 机制详解

### 什么是 Skill

Skill（技能）是 OpenCode 的可复用指令机制，允许你通过 `SKILL.md` 文件定义可复用的行为。技能通过原生的 `skill` 工具按需加载——代理可以查看可用技能，并在需要时加载完整内容。

### Skill 的放置位置

OpenCode 会搜索以下位置来发现技能：

- 项目配置：`.opencode/skills/<name>/SKILL.md`
- 全局配置：`~/.config/opencode/skills/<name>/SKILL.md`
- 项目 Claude 兼容：`.claude/skills/<name>/SKILL.md`
- 全局 Claude 兼容：`~/.claude/skills/<name>/SKILL.md`

### Skill 的结构

每个 `SKILL.md` 必须以 YAML frontmatter 开头，包含以下必填字段：

- `name`：技能名称（1-64个字符，仅小写字母和数字，可用单个连字符分隔）
- `description`：技能描述（1-1024个字符）

示例：

```markdown
---
name: web-search
description: 用于网络搜索的技能
license: MIT
---

## What I do
- 根据用户需求进行网络搜索
- 获取并总结搜索结果

## When to use me
当用户需要查询实时信息或最新资讯时使用此技能
```

### 如何使用 Skill

在配置阶段，可以通过设置来决定代理需要预先加载何种技能。所谓加载，其实就是把SKILL.md中的内容加载到上下文。

```
skill({ name: "web-search" })
```

另外，大部分技能也描述了什么时候会调用，因此如果对话的提示词命中了工具被调用的状况，也会自动加载技能（除非定义时被强制deny了）

## 实际使用

在我这次的工作中，大致的步骤可以拆成3部分：

- 信息检索和收集
- 话题选择
- 信息验证

然后再来一个主agent用来控制前述工作流。

## 遇到的问题

### 来源链接失真

这个问题比我预想的还严重。

我让Agent先通过Agent-Browser来检索信息，然后根据检索到的信息来生成内容，这些部分其实完成的很不错，输出的内容并不会太失真。

可是唯独保存信息来源链接这个看上去很简单的事情，初开发时用的Minimax-M2.5就是做不到，保存下来的连接10个9个都是错的。即使将Temprature设置到最低，附加更多的提示词，也无法保存真实的链接。

和对信息的时候，情况更是让人哭笑不得，它保存链接的时候，极大概率会修改链接中单词的单复数，和一些代号字符串。比如:

- video > videos
- a-s001 > b-a273

最后折腾了半天，换成GLM就没有这个问题了...

### Agent动作不稳定

这个问题涉及生成过程的方方面面。

对于主控Agent，当它负责调度多个subgent时，其主控稳定性并没有那么灵活，比如我的工作流有3个步骤，有的时候我想只执行其中的一部分，但是我输入的指令只是对话中Prompt的一小段问题，模型最终很可能还是根据更长的预定义Propmt把整个流程运行完了。

对于Subgent，一旦离开了之前的上下文，或没有明确指定实现方式，它每一的实现方式，可能会完全不相同。比如这一次写js，下一次写python，再下一次直接Shell执行。整个过程像极了早期文生图模型的反复抽卡...

## 获得的经验

- 不同模型擅长点不同，需要配合使用
    + GLM：指令执行准确，但是输出非常慢
    + Minimax：非Coding的指令执行不稳定，但写代码不错，输出很快
    + kimi：指令执行写代码都适中，但是会莫名其妙死循环，不适合用在控制Agent上

- 确定性高的步骤，使用Agent的效率不如直接生成代码和Skill，让后让Agent来直接调用

- 工作流还是有其意义，相比基于强Agent的多Agent体系，其稳定性更好，代价就是牺牲一些灵活性，需要根据具体任务来选择。