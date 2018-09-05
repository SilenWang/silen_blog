---
title: SGE系统使用笔记
categories: Bioinfomatic
date: 2018-08-30 09:36:38
tags:
---

在新公司使用SGE的方式有点不太一样, 做个笔记

<!-- more -->

# 脚本编写

脚本写完后, 并不在投递时进行参数指定, 而是写到被投递的脚本中, 然后直接`qsub shell.sh`

```bash
#！/bin/bash
#$ -S /bin/bash          //表明此脚本为bash
#$ -V                    //传递当前命令的所有环境变量
#$ -cwd                  //将当前路径设置为工作路径
#$ -N WorkName           //任务名
#$ -o WorkName.log       //任务输出日志文件名
#$ -j y                  //指定任务的标准错误流是否合并到标准输出流y[es] n[o]

shell script
```

投递含有上述内容本质上等价于:

```bash
qsub -S -V -cwd \
    -N WorkName \
    -o WorkName.log \
    -j y \
    shell script
```

# qlogin的使用

在进行测试时, 可使用`qlogin`登陆到计算节点, 然后在计算节点环境下进行的操作与写好脚本`qsub`投递到集群上效果一致