---
title: 在Jupyter中注册新的Python或R内核
categories: Script
date: 2023-01-17 12:53:40
tags: ['jupyter']
---

进行测试时经常要新建一个conda环境, 然后将环境装上jupyter内核后用notebook进行测试. 每次注册新内核都要重查, 所以还是做个记录好了...
<!-- 摘要部分 -->
<!-- more -->

## Python内核准备和注册

- 进入要注册的环境中后, 再进行注册

```bash
conda activate YOUR_ENV
conda install jupyter ipykernel
python -m ipykernel install --user --name ENV_NAME
```

## R内核准备和注册

- 进入要注册的环境中后, 再进行注册

```bash
conda activate YOUR_ENV
conda install jupyter r-irkernel
```

```r
IRkernel::installspec(name = 'REG_NAME', displayname = 'ENV_NAME')
```