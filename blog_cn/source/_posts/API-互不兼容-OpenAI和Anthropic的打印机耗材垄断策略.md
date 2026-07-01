---
title: API 互不兼容？OpenAI 和 Anthropic 的"打印机耗材垄断"策略
date: 2026-07-01 22:00:00
tags:
- AI
- API
- OpenAI
- Anthropic
- 吐槽
categories: Others
---

最近在做一些 AI 应用集成的工作，需要在不同模型之间切换。然后我就发现了一个让人非常抓狂的事实——OpenAI 和 Anthropic 的 API 几乎是**完全不兼容**的。

<!-- more -->

## 从一次简单的换模型说起

假设你已经在用 OpenAI 的 GPT-4o，现在想试试 Claude 3.5 Sonnet。按照直觉，换个 endpoint、改个 model 名字就行了吧？天真。

```python
# OpenAI 的调用方式
import openai
client = openai.OpenAI(api_key="sk-xxx")
response = client.chat.completions.create(
    model="gpt-4o",
    messages=[{"role": "user", "content": "Hello"}]
)
```

而 Anthropic 呢？

```python
# Anthropic 的调用方式
import anthropic
client = anthropic.Anthropic(api_key="sk-ant-xxx")
response = client.messages.create(
    model="claude-3-5-sonnet-20241022",
    max_tokens=1024,
    messages=[{"role": "user", "content": "Hello"}]
)
```

一眼看去好像差别不大，但实际用起来处处是坑：

- **客户端 SDK 不同**：一个是 `openai` 包，一个是 `anthropic` 包，没法混用
- **参数命名不同**：`max_tokens` 的位置和写法不一样
- **返回格式不同**：OpenAI 返回流式用 `for chunk in response`，Anthropic 用 `for event in response` 还得判断 event type
- **stream 实现细节不同**：连处理流式响应都得重写
- **Tool/function calling 协议不同**：这是最要命的——如果你的应用用了 function calling，切换模型基本等于重写整个调用层
- **系统提示词格式不同**：OpenAI 是 `system` role，Anthropic 有自己的 system prompt 参数
- **图像输入格式不同**：连传张图都得用不同格式

这种差异逼着开发者要么写一堆适配代码，要么用第三方库（如 LiteLLM、OpenRouter）做一层封装。

## 这让我想起了打印机

这场景太熟悉了。当年买打印机，机器本身可以卖得很便宜，但墨盒/硒鼓一定是专用的。HP 的打印机只能用 HP 的墨盒，Canon 的只能用 Canon 的。第三方兼容墨盒要么被固件封锁，要么被"官方认证芯片"挡住。

换个角度想想：
- **OpenAI = HP**，家大业大，生态完善，大家都用它的 API 格式
- **Anthropic / Claude = Canon**，后来者，性能不错，但偏要自己搞一套协议
- **你 = 被锁定的用户**，一旦深度集成了一家的 API，迁移成本高得吓人

## API 锁定的三重陷阱

OpenAI 和 Anthropic 的 API 不兼容不完全是"技术债"或"设计理念不同"，背后有很清晰的商业逻辑：

1. **切换成本**：你在一家 API 上写了大量代码、做了深度集成，就很难说换就换。即使另一家模型更好、更便宜，心理和工程上的迁移成本都在那里。
2. **生态锁定**：围绕 OpenAI API 已经长出了庞大的生态——各类工具、监控、代理、缓存层。换 API 意味着这些东西都得跟着改。
3. **数据壁垒**：虽然 API 不直接锁数据，但当你深度依赖某个平台的特性（比如 OpenAI 的 structured output、Anthropic 的 extended thinking），就被绑得更死了。

这不就是数字时代的打印机耗材吗？机器（API 标准）本身可以开放，但耗材（实现细节、生态工具、专属特性）是专用的。

## 好在有一些破局者

当然，社区也不是吃素的。已经有不少项目在尝试打破这种锁定：

- **LiteLLM**：提供统一接口，覆盖 100+ 模型供应商
- **OpenRouter**：一个 API 调用多个模型
- **Vercel AI SDK**：框架级的抽象层
- **LangChain**：虽然被吐槽过度封装，但这种场景下确实有用

但这些方案本质上都是在 OpenAI 的 API 风格上做适配——Anthropic 的原始 API 并没有因此变成"标准"。

## 太阳底下无新事

从打印机到智能手机充电口，从游戏主机独占到云服务厂商的锁定策略，商业公司通过制造不兼容来绑架用户，这招已经用了上百年。

AI 时代也不例外。OpenAI 和 Anthropic 的 API 互不兼容，说白了就是——**当年打印机厂商那套"机器随便买，耗材必须用我的"策略，换个马甲又回来了**。

作为开发者，我们能做的就是：
1. 在设计系统时预留抽象层，不要把代码和某家 API 绑死
2. 关注标准化进展，比如有没有统一的 API 标准出现
3. 支持那些致力于开放生态的项目

不过说实话，指望商业公司主动放弃锁定策略，就跟指望打印机厂商免费送墨盒一样——想得美。

太阳底下无新事。
