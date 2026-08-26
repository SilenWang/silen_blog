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

DeepSeek 最近正式涨价了。作为一个上季度 Token 消耗直奔上亿的家伙（详见 {% post_link 三个月Token消耗暴涨13倍-我的AI工作流演变记录 [《三个月Token消耗暴涨13倍》] %}），这消息本来我没太在意的，毕竟 DeepSeek 本来就非常便宜，涨价后的价格，其实还是比V4最开始宣布时的正式价格低的。但是我还是低估了我的用量，涨价后一天的实测下，这价格确实个人难以承受，不想想办法要付费上班了。

<!-- more -->

## Token Plan 的"最甜档位"悖论

为了省钱，我开始研究各家卖 Token Plan（包月/包量套餐）的商家，结果发现一个很有意思的现象：

> 几乎每家都有一个最值的档位——按 Token 单价算，最为划算；而是更贵的套餐，单价其实反而会稍微上升。

这大概是商家故意设计的：低价档位用来引流，买的人多、单价亏一点无所谓；大额套餐则是给"懒得折腾"的用户准备的，单价自然要高一些。换句话说，**无脑买最贵的套餐，不一定合算**。

所以正确的省钱姿势就很清晰了：

1. 多家商家，每家都买它的"最值档位"
2. 注册多个账户，把总量凑够
3. 反正 DeepSeek、Kimi、GLM、MiniMax 这些模型大家都有，谁家便宜用谁家

但这么做有个巨大的麻烦：**Key 太多了**。假设你有 5 家 × 2 个账户 = 10 个 Key，而你的工具链里还有 Claude Code、Codex、opencode……每个工具都要配 Key、换 Key、记住哪家还剩多少额度。光是管理这些 Key 就能把人逼疯。

## Omniroute：一个 Key 管所有

这时候就需要 Omniroute 这种项目了。

[Omniroute](https://github.com/diegosouzapw/OmniRoute) 是一个开源的 AI 网关（MIT 协议，GitHub 上已经 5 万多 star），干的事一句话总结：**把几百家模型供应商聚合成一个统一入口，客户端只需要一个 Key**。

它的特性它自己的官网有详细的对比，简单总结就是商家支持广泛，路由策略丰富。

用起来的大概流程：先部署一个 Omniroute（我用的 Docker），把各个商家的 Key 填进控制台，然后用Agent调用MCP来配置策略
配置好路由策略（比如"优先用最便宜的，额度不够自动回退"），再把所有 AI 工具的 `base_url` 指过去。之后想切商家？控制台里点一下就行，客户端完全无感。

Omniroute 还做了不少易用性改进，比如它带一个 **MCP server**（暴露了 100 多个工具），可以让 AI 自己帮你配置——切换组合、设置路由策略、查询额度，说句话就行。对于我这种懒人来说，确实省了不少事。

## 实际部署中踩到的坑

不过实际用下来，Omniroute 不是一点问题没有。

### 坑一：搜索拦截功能有问题，Codex 和 Claude Code 搜不了网

Omniroute 有个"搜索拦截"功能：把客户端自带的搜索工具（比如 Claude Code 的 `web_search_20250305`）拦截下来，转发到它自己的 `/v1/search` 网关，用你配置的搜索供应商（Exa、Serper 这些）去执行。

这个功能其实很好，因为codex和claude code都会需要服务端提供搜索结果，而国外很多服务商只提供token，不提供工具那些的。但这功能在目前还不能实际使用。Claude Code 发的是带版本号的 `web_search_20250305`，而 Omniroute 的拦截匹配器只匹配裸的 `web_search`——版本号对不上，拦截根本不生效，上游模型返回一个空的搜索结果，Claude Code 拿到手就是"搜了个寂寞"。Codex 那边也类似，它的独立搜索接口对自定义供应商不支持。

这个坑影响说大不大，毕竟用别的方式，技能或者mcp让agent能获得结果就行了。但是说小也不小，本来是打算omniroute解决一切问题的，现在只能等修复了。相关 issue 目前还还挂着（[#9772](https://github.com/diegosouzapw/OmniRoute/issues/9772)、[#9725](https://github.com/diegosouzapw/OmniRoute/issues/9725)、[#8674](https://github.com/diegosouzapw/OmniRoute/issues/8674)）。

### 坑二：太新的商家同步不了额度

比如我用的 Commandcode Go，推理完全正常，但额度在控制台的 Provider Quota 里就是不显示。更麻烦的是，额度同步不上可能会导致"额度感知的自动回退"失效：账户额度耗尽了不会被正确标记，多账户轮换机制就乱了，甚至可能出现所有账户都被打上"额度耗尽"标签、重置后也不自动恢复的情况（[#9603](https://github.com/diegosouzapw/OmniRoute/issues/9603)）。

也就是说，**越是新出的、越便宜的商家，越容易遇到额度同步问题**，而这恰恰是我们最想用的那些。

### 坑三：默认不允许重型任务并发，用着用着 Agent 突然就断了

前两个坑好歹能忍，这个坑是真能气死人：**在默认配置下，Omniroute 并不允许重型任务并发（OMNIROUTE_CHAT_MAX_HEAVY_IN_FLIGHT=1），并行跑超过一个任务 Agent 会突然断掉。**

现在的 Agent 基本都会一个大任务会同时拆出好几个子任务，并行地打一堆 LLM 请求——比如 opencode 的子 Agent、Claude Code 的并行工具调用。这一下并发就上去了。

## 小结

虽然不圆满，但现在这个方案少说也能解决80%时候的需求了，希望前面的问题能早日更新吧。。。

最后照例附上项目地址：[github.com/diegosouzapw/OmniRoute](https://github.com/diegosouzapw/OmniRoute)，有需要的朋友可以自行研究。
