---
title: I always stumble on whatever I use – several issues with Devcontainer
categories: Others
date: 2025-12-19 22:27:38
tags: ['devcontainer', 'devpod', 'docker', 'git']
---

I just wanted to set up a devcontainer environment to maintain the company website more efficiently, but I never expected to encounter three pitfalls in a single task...

<!-- more -->

## Background

[Devcontainer](https://code.visualstudio.com/docs/devcontainers/containers) is a feature introduced by Visual Studio Code that allows us to quickly set up a development environment matching the project using Docker containers. This "code-as-environment" approach can greatly unify the team collaboration experience and reduce the "it works on my machine" embarrassment. The devpod I wrote earlier is an open‑source solution based on this tech.

However, the ideal is plump, while reality always trips you up on the details. While configuring a Devcontainer for my company’s official website project recently, I ran into three unexpected problems in a row, each costing a fair amount of time to solve. Below I’ll document these issues and their solutions for future reference.

## Pitfall 1: Missing Git subtree command

### Problem description
The base Devcontainer images provided by Microsoft (Debian and Ubuntu) both come with Git, but the version is quite basic. It does **not** include the `subtree` subcommand… Well, this was the first time I learned that Git has so many features that some relatively new commands are made optional…

But that’s not the real problem. The Git packages in Ubuntu and Debian repositories should be fully functional—outside the Devcontainer I’ve never seen a “subtree command not found” message. The issue appeared even after trying `apt-get install git-man` or `git-all`; the `subtree` functionality still wasn’t available. Why? Actually it’s simple: the Git that comes inside the Devcontainer is **not** from the distribution repositories; it’s located in `/usr/local/bin`, while the system‑installed Git would be placed in `/usr/bin`. Since `/usr/local/bin` has higher priority in the `PATH`, the Git we install via `apt` never gets used, and `subtree` remains missing.

### Solution
Given the above, the fix is straightforward: install the distribution’s Git via `apt` and then remove the one that came with the image:

```dockerfile
FROM mcr.microsoft.com/devcontainers/base:ubuntu

RUN apt-get update && \
   add-apt-repository -y ppa:git-core/ppa && \
   apt-get update && \
   apt-get install -y git vim && \
   # remember to remove the original git from the container
   rm -r /usr/local/bin/git* && \
   apt-get autoremove -y && \
   apt-get clean -y
```

## Pitfall 2: Automatic line‑ending conversion on Windows

### Problem description
When you run a Devcontainer on a Windows host using Docker Desktop, Git will automatically convert line endings from LF to CRLF according to the default setting of `core.autocrlf` (usually `true`). This is a common behavior in Windows development environments, but inside the container the file system is Linux‑style, so Git performs the conversion during checkout, causing almost every file to be marked as “modified”.

The concrete symptom is: executing `git status` inside the container shows a large number of files as changed (even though you did nothing). This noise is not only confusing but may also affect subsequent commit operations.

### Solution
Well… I didn’t actually solve this problem, because tweaking Git configuration doesn’t stop VS Code from performing the conversion… This issue only occurs when using Docker Desktop on Windows. Anyway, compiling a website with Docker Desktop can even freeze the whole desktop, so just don’t use it… Stick to using a linux server as an SSH Provider.

## Pitfall 3: Devcontainer Desktop and remote SSH Provider disconnections

### Problem description
When trying to connect to a remote SSH Provider through Devpod Desktop, the connection turns out to be extremely unstable—it often drops within a few minutes, and the reconnection process is painfully slow. It disconnects again shortly after re‑establishing. Since Devpod Desktop doesn’t provide any runtime logs, I have no clue what’s wrong…

### Solution
Simple and brutal: use Devpod Desktop only for configuring Providers and workspaces, and let the DevPod CLI handle the actual connection.

## Summary
Good news: all problems were eventually solved.  
Bad news: I lost another seven or eight hours of sleep this week. May I not drop dead in 2025…
