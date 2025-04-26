---
title: Using aider for programming assistance
categories: Script
date: 2025-04-25 00:42:46
tags: ['AI', 'aider']
---

This year, DeepSeek has made LLMs popular again. Although after ChatGPT, Copilot and Cursor have been advertising everywhere, I never actually tried them... So this year, I... tried [Cline](https://cline.bot/), which is indeed good, but also matches my preconception - not that convenient. Especially since I need something that can be used directly in the command line and quickly integrated into scripts to batch complete some simple coding initialization tasks. That's how I found [aider](https://aider.chat/)...

<!-- more -->

## Installation

The installation is quite simple since it's a Python project, just use pip:

```bash
python -m pip install aider-install
aider-install
```

## Usage

aider is a CUI application, just enter `aider` to start. Aider supports a series of chat commands starting with `/` that help you more efficiently manage code editing tasks, switch modes, and interact with the codebase. Here are some commonly used commands and their functions:

- File management commands
    + `/add`: Add files to the chat session so Aider can edit or review them in detail.
    + `/drop`: Remove files from the chat session to free up context space.
    + `/ls`: List all known files and indicate which ones are included in the chat session.
- Editing and code operation commands
    + `/code`: Request code changes, or switch to code mode if no prompt is provided.
    + `/diff`: Show the diff of changes since the last message.
    + `/commit`: Commit changes made to the codebase outside the chat (optional commit message).
    + `/undo`: Undo Aider's last Git commit.
- Mode switching commands
    + `/chat-mode`: Switch to new chat mode.
    + `/architect`: Enter architect/edit mode, using two different models for code design or editing.
    + `/ask`: Ask questions about the codebase without editing any files.
- Other utility commands
    + `/clear`: Clear chat history.
    + `/copy`: Copy the assistant's last message to clipboard.
    + `/model`: Switch to a new language model.
    + `/help`: Get help about using Aider.
    + `/exit` or `/quit`: Exit the Aider application.

## Embedding in shell scripts

As a CUI program, aider can also be run directly from the command line by specifying parameters and writing the prompt directly in the command line (this is quite similar to `ssh`).

### Use Case 1 - Batch blog translation

Since it can be run from the command line, aider can be conveniently embedded in shell scripts to batch complete some not-too-heavy tasks, like translating all past blogs...

Although when ChatGPT first became popular, I thought about using it for translation, but after actual operation I found there was still significant manual work. Because I already have hundreds of blog files... Each time entering the prompt, then copying content in and out and proofreading was still tedious and boring work.

Using aider, I can at least be a bit lazier.

```bash
set -e 

DIR_A=blog_cn/source/_posts
DIR_B=blog_en/source/_posts

for file in `ls $DIR_A/*`; do
    # Get filename (without path)
    filename=$(basename "$file")

    # Check if file exists in directory B
    if [[ -e "${DIR_B}/${filename}" ]]; then
        echo "File ${filename} exists in ${DIR_B}, skipping"
    else
        echo "File ${filename} doesn't exist in ${DIR_B}, starting aider translation..."
        echo "Starting document proofreading"
        aider \
            --no-show-model-warnings --yes --no-auto-commits \
            --message "I have a blog post here, please translate it to English and save it to ${DIR_B}/${filename}" \
            $file
    fi
```

### Use Case 2 - Batch API test generation

Similar to blog generation, theoretically test code for interfaces can be quickly written based on OpenAPI or Swagger documentation, but if there are many interfaces, writing them manually is still troublesome. In such cases, having aider read the code and pre-generate a test code framework can save a lot of effort.

```bash
for group in group1 group2 group3;do
    # Get filename (without path)
    pytest=test/${group}.py
    touch $pytest
    aider \
        --no-show-model-warnings --yes --no-auto-commits \
        --message "Here is an API's Swagger description file ${group}.swagger.json, please write pytest test scripts based on its content and save the test scripts to ${pytest}." \
        data/${group}.swagger.json $pytest
done
```
