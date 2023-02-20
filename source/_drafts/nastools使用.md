---
title: nastools使用
categories: Daily
date: 2023-02-04 00:46:30
tags: ['nastools', 'docker', 'emby']
---


<!-- 摘要部分 -->
<!-- more -->

- 注册TMDB会员, 获取API_Key(114c56e983b64425b2a51b46ccaf64b5)
- 关注BotFater获取Token(6008501130:AAEEzHwLwXPwkwvpH-vkvPBZeCsR2mTf5_o)
- 准备目录结构
- 拉取docker, 编写`docker-compose.yaml`
- 获取`emby`镜像 `zishuo/embyserver`
  + 启动镜像, 获取API KEY, 设置WEBHOOK 
  + 将KEY填到Nastools的配置中
- 获取下载软件镜像 `linuxserver/qbittorrent`
  + 启动, 更改密码
  + 在Nastools进行设置 