---
title: 使用aider辅助编程工作
categories: Script
date: 2025-04-25 00:42:46
tags: ['AI编程', '命令行工具', '自动化', '批量处理', 'aider', 'AI-assisted programming', 'command-line tool', 'automation']
---

今年的deepseek着实让LLM又火了一把，虽然在chatGPT之后，copilot和cursor就一直在到处打广告，但是我其实一直没试过... 于是在今年，我... 试了试[Cline](https://cline.bot/)，确实好用，但是也符合我的固有印象，没有那么方便，尤其是，我更需要一个能直接在命令行使用并能够快速迁入脚本的选项，以批量完成一些简单的编码初始化业务，于是我找到了[aider](https://aider.chat/)...

<!-- more -->

## 安装

软件安装比较简单，因为是一个python项目，使用pip就可以安装：

```bash
python -m pip install aider-install
aider-install
```

## 使用

aider是一个CUI应用，输入aider就能进入。Aider 支持一系列以 `/` 开头的聊天命令，这些命令可以帮助你更高效地管理代码编辑任务、切换模式以及与代码库进行交互。以下是一些常用的命令及其功能：
- 文件管理命令
    + `/add`：将文件添加到聊天会话中，以便 Aider 可以编辑或详细审查这些文件。
    + `/drop`：从聊天会话中移除文件，释放上下文空间。
    + `/ls`：列出所有已知文件，并指示哪些文件已包含在聊天会话中。
- 编辑与代码操作命令
    + `/code`：请求对代码进行更改，如果没有提供提示信息，则切换到代码模式。
    + `/diff`：显示自上次消息以来的更改差异。
    + `/commit`：提交在聊天之外对代码库所做的更改（可选提交信息）。
    + `/undo`：撤销 Aider 执行的最后一次 Git 提交。
- 模式切换命令
    + `/chat-mode`：切换到新的聊天模式。
    + `/architect`：进入架构师/编辑模式，使用两种不同的模型进行代码设计或编辑。
    + `/ask`：询问有关代码库的问题，而不编辑任何文件。
- 其他实用命令
    + `/clear`：清除聊天历史记录。
    + `/copy`：将上次助手的消息复制到剪贴板。
    + `/model`：切换到新的语言模型。
    + `/help`：获取有关 Aider 使用的帮助信息。
    + `/exit` 或 `/quit`：退出 Aider 应用程序。

## 嵌入shell脚本

作为CUI程序，aider 也有直接从命名行运行的方式，通过指定参数，然后直接将 prompt 写在命名行中，就可以（这方式还挺类似`ssh`）。

### 用例1--批量翻译博客

因为可以通过命令行运行，aider也就可以方便的潜入到shell中，批量完成一些工作量不太大的问题，比如将过过去的博客都翻译一遍...

虽然在chatGPT刚火的时候，我就想过用它来完成翻译，但是实际操作后才发现，手动工作量仍然不小。因为我的博客已经有上百个文件了... 每次输入prompt,然后把内容贴近贴出以及进行校对，也还是很繁琐和无趣的工作。

使用aider，至少能再多偷些懒。

```bash
set -e 

DIR_A=blog_cn/source/_posts
DIR_B=blog_en/source/_posts

for file in `ls $DIR_A/*`; do
    # 获取文件名（去掉路径部分）
    filename=$(basename "$file")

    # 检查文件是否存在于目录B
    if [[ -e "${DIR_B}/${filename}" ]]; then
        echo "文件 ${filename} 存在于目录${DIR_B}，跳过"
    else
        echo "文件 ${filename} 不存在于目录${DIR_B}，启动aider翻译..."
        echo "开始文稿校对"
        aider \
            --no-show-model-warnings --yes --no-auto-commits \
            --message "我这里有一篇blog，请把它翻译成英文，然后保存到 ${DIR_B}/${filename} 文件中" \
            $file
    fi
```

### 用例2--批量生成api测试

与生成博客类似，根据openapi或者swagger文档理论上可以很快的编写出接口的测试代码，但是如果接口一多，实际写起来还是很麻烦的。这种时候，让aider来读取代码，并预先生成一个测试代码的框架，可以省下不少事情。

```bash
for group in group1 group2 group3;do
    # 获取文件名（去掉路径部分）
    pytest=test/${group}.py
    touch $pytest
    aider \
        --no-show-model-warnings --yes --no-auto-commits \
        --message "这里有一份api的swagger描述文件 ${group}.swagger.json ，请根据它的内容，为这个api编写pytest测试脚本，并将写好的测试脚本放到 ${pytest} 文件中。" \
        data/${group}.swagger.json $pytest
done
```
