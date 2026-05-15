---
title: Installing and Using GitHub CLI on Ubuntu
categories: Tools
date: 2026-05-10 12:00:00
tags: ['GitHub CLI', 'gh', 'Ubuntu', 'Authentication', 'Command Line', 'git', 'issue']
---

[GitHub CLI](https://cli.github.com/) (`gh`) is GitHub's official command-line tool that lets you interact with GitHub without leaving your terminal. Even better, once authenticated, `git` commands automatically reuse these credentials — no SSH configuration needed.

<!-- TLD; DR -->
<!-- more -->

## Installation Methods

### Method 1: Install via apt (Recommended)

GitHub CLI provides an official apt repository:

```bash
# Install dependencies
sudo apt install curl

# Add GitHub CLI's GPG key and repository
curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | sudo dd of=/usr/share/keyrings/githubcli-archive-keyring.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" | sudo tee /etc/apt/sources.list.d/github-cli.list > /dev/null

# Update and install
sudo apt update
sudo apt install gh
```

Verify the installation:
```bash
gh --version
```

### Method 2: Install via snap

If you prefer snap packages:

```bash
sudo snap install gh
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

```bash
# Clone a private repo (no extra config needed)
git clone https://github.com/yourname/private-repo.git

# Push code (reuses gh's credentials automatically)
git push origin main
```

## Common GitHub CLI Commands

### Working with Issues

```bash
# List issues in the current repository
gh issue list

# View issue details
gh issue view 123

# Create a new issue
gh issue create --title "bug report" --body "describe the bug"

# Close an issue
gh issue close 123
```

### Working with Pull Requests

```bash
# List PRs
gh pr list

# Create a PR
gh pr create --title "fix bug" --body "fix description"

# View PR details
gh pr view 123

# Merge a PR
gh pr merge 123
```

### Repository Operations

```bash
# View repository info
gh repo view

# Open repository in the browser
gh repo view --web

# Fork a repository
gh repo fork
```

## Summary

GitHub CLI makes working with GitHub from the terminal incredibly smooth. A single `gh auth login` handles all authentication, and after that, both `gh` commands and regular `git` operations work without credential prompts. For any developer who frequently interacts with GitHub, this tool is well worth installing.
