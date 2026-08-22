---
title: Omniroute：DeepSeek涨价后，我用一个Key管起了所有Token Plan
categories: AI
date: 2026-08-22 22:00:00
tags:
  - Omniroute
  - DeepSeek
  - Token
  - AI
  - API
  - 省钱
---

DeepSeek 最近又涨价了。作为一个上个月 Token 消耗直奔上亿的家伙（详见 {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [《三个月Token消耗暴涨13倍》] %}），这消息直接把我干沉默了：用量省不下来，那就只能在单价上做文章了。

<!-- more -->

## 涨价的 DeepSeek，和上亿的 Token 消耗

先说背景。之前 DeepSeek V4-Pro 刚发布的时候，官宣永久降价到原价四分之一，我当时还专门写了篇文章夸它良心（{% post_link DeepSeek-永久降价-大模型价格战的新玩法 [《DeepSeek 终于还是调价了... 等等... 降价？》] %}）。结果才过了几个月，价格又涨回去了……

而我的消耗量早就不是当初那个量级了。上篇文章统计过，三个月 Token 消耗从 38M 暴涨到 500M/月，最近更是稳定在一个月上亿的量级。用量上去了、单价也上去了，这账没法看。

## Token Plan 的"最甜档位"悖论

为了省钱，我开始研究各家卖 Token Plan（包月/包量套餐）的商家，结果发现一个很有意思的现象：

> 几乎每家都有一个"最甜"的档位——按 Token 单价算，它往往是最便宜的；反而是越贵的套餐，单价越贵。

这大概是商家故意设计的：低价档位用来引流，买的人多、单价亏一点无所谓；大额套餐则是给"懒得折腾"的用户准备的，单价自然要高一些。换句话说，**无脑买最贵的套餐，反而是最亏的**。

所以正确的省钱姿势就很清晰了：

1. 多家商家，每家都买它的"最甜档位"
2. 注册多个账户，把总量凑够
3. 反正 DeepSeek、Kimi、GLM、MiniMax 这些家的模型能力大差不差，谁家便宜用谁家

但这么做有个巨大的麻烦：**Key 太多了**。假设你有 5 家 × 2 个账户 = 10 个 Key，而你的工具链里还有 Claude Code、Codex、opencode……每个工具都要配 Key、换 Key、记住哪家还剩多少额度。光是管理这些 Key 就能把人逼疯。

## Omniroute：一个 Key 管所有

这时候就需要 Omniroute 这种项目了。

[Omniroute](https://github.com/diegosouzapw/OmniRoute) 是一个开源的 AI 网关（MIT 协议，GitHub 上已经 5 万多 star），干的事一句话总结：**把几百家模型供应商聚合成一个统一入口，客户端只需要一个 Key**。

几个关键特性：

- **支持的商家极多**：340+ 家供应商（其中 90+ 家还有免费额度）、1200+ 模型，Kimi、Claude、GPT、Gemini、GLM、DeepSeek、MiniMax 这些主流的全覆盖
- **聚合策略丰富**：内置 19 种路由策略，支持四层自动降级（订阅 → API → 便宜 → 免费）、额度感知的自动回退、熔断、Key 冷却
- **使用端零改动**：把 Claude Code / Codex / Cursor 这些工具的 `base_url` 指向 `http://localhost:20128/v1` 就行，工具侧只认一个 Key，剩下的全在 Omniroute 控制台里切换
- **还能省 Token**：内置 RTK + Caveman 压缩引擎，号称能省 15%~95% 的 Token
- **部署方式随意**：npm 全局安装、Docker、Electron 桌面版、PWA 都行

用起来的大概流程：先部署一个 Omniroute（我用的 Docker），把各个商家的 Key 填进控制台，配置好路由策略（比如"优先用最便宜的，额度不够自动回退"），再把所有 AI 工具的 `base_url` 指过去。之后想切商家？控制台里点一下就行，客户端完全无感。

Omniroute 还做了不少易用性改进，比如它带一个 **MCP server**（暴露了 100 多个工具），可以让 AI 自己帮你配置——切换组合、设置路由策略、查询额度，说句话就行。对于我这种懒人来说，确实省了不少事。

## 实际部署中踩到的坑

不过实际用下来，Omniroute 也远没到"开箱即用"的程度，至少我踩到了两个比较闹心的坑。

### 坑一：搜索拦截功能有问题，Codex 和 Claude Code 搜不了网

Omniroute 有个"搜索拦截"功能：把客户端自带的搜索工具（比如 Claude Code 的 `web_search_20250305`）拦截下来，转发到它自己的 `/v1/search` 网关，用你配置的搜索供应商（Exa、Serper 这些）去执行。

想法很好，但实际有 bug。Claude Code 发的是带版本号的 `web_search_20250305`，而 Omniroute 的拦截匹配器只匹配裸的 `web_search`——版本号对不上，拦截根本不生效，上游模型返回一个空的搜索结果，Claude Code 拿到手就是"搜了个寂寞"。Codex 那边也类似，它的独立搜索接口对自定义供应商不支持。

这个坑影响挺大：本来聚合的初衷是"一个 Key 管所有"，结果搜索功能直接废了，相当于把 AI 工具的联网能力砍掉一半。相关 issue 还挂着（[#9772](https://github.com/diegosouzapw/OmniRoute/issues/9772)、[#9725](https://github.com/diegosouzapw/OmniRoute/issues/9725)、[#8674](https://github.com/diegosouzapw/OmniRoute/issues/8674)），暂时只能等修复。

### 坑二：太新的商家同步不了额度

Omniroute 的"额度同步"是白名单制的——每个商家的额度查询接口都要单独注册，新商家经常漏掉。

比如我用的 Qwen Cloud Token Plan，推理完全正常，但额度在控制台的 Provider Quota 里就是不显示。更麻烦的是，额度同步不上会导致"额度感知的自动回退"失效：账户额度耗尽了不会被正确标记，多账户轮换机制就乱了，甚至可能出现所有账户都被打上"额度耗尽"标签、重置后也不自动恢复的情况（[#9603](https://github.com/diegosouzapw/OmniRoute/issues/9603)）。

也就是说，**越是新出的、越便宜的商家，越容易遇到额度同步问题**，而这恰恰是我们最想用的那些。

## 总结

用了一段时间，我的感受是：

- Omniroute 解决了核心问题——**一个 Key 管所有账户**，聚合策略也确实好用，省心不少
- 但它毕竟还年轻，**搜索拦截、额度同步这类"外围功能"还有不少坑**，而且越是新商家越容易踩
- 如果只是把它当"多 Key 聚合器"用，体验很好；但如果依赖它的搜索网关、或者指望所有商家都能同步额度，建议先看看对应 issue 的进展

最后照例附上项目地址：[github.com/diegosouzapw/OmniRoute](https://github.com/diegosouzapw/OmniRoute)，有需要的朋友可以自行研究。
