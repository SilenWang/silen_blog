---
title: Installing and Using GitHub CLI on Ubuntu
categories: Tools
date: 2026-05-10 12:00:00
tags: ['GitHub CLI', 'gh', 'Ubuntu', 'Authentication', 'Command Line', 'git', 'issue']
---

Every time I get a new machine, I used to configure SSH keys for GitHub to clone and push code. But I never bothered to keep track of the keys, so I ended up recreating them repeatedly. It wasn't until I ran into `multica` and browsed the issue tracker that I realized GitHub provides GitHub CLI (`gh`) for quick authentication...

<!-- more -->

## Installation Methods

### Method 1: Install via apt (Recommended)

Debian 13 and newer versions of Ubuntu have included `gh` in their official repositories, so you can install it directly:

```bash
sudo apt update
sudo apt install gh
```

## Authentication

After installation, the first step is to log into your GitHub account:

```bash
gh auth login
```

This starts an interactive flow:

1. Choose **GitHub.com** as the account type
2. Choose **HTTPS** as the protocol
3. Choose **Login with a web browser** and copy the one-time code displayed
4. Paste the code in your browser and authorize
5. The terminal will confirm successful authentication

For non-interactive authentication using a Token:

```bash
gh auth login --with-token < my_token.txt
```

## How Git Benefits from `gh` Authentication

After `gh auth login`, `git` automatically uses `gh`'s credential manager:

- `git push`, `git pull`, and other operations no longer ask for username/password for HTTPS remotes
- No need to manually configure SSH keys for code operations
- All authentication is centrally managed through `gh`
