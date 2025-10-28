---
title: Setting Up Self-Hosted Codespaces with DevPod
categories: Coding
date: 2025-10-28 05:49:54
tags: [DevPod, DevOps, Development Environment, Docker, Podman]
---

DevPod is an open-source development environment management tool that allows you to create development environments similar to GitHub Codespaces on any Kubernetes cluster or Docker host. This article will introduce how to use the DevPod CLI to create workspaces and provide a detailed explanation of writing DevContainer configuration files.

<!-- more -->

## Introduction to DevPod

DevPod is an open-source tool developed by Loft Labs that enables developers to create reproducible, disposable development environments on any infrastructure. Similar to GitHub Codespaces, but DevPod is self-hosted and can run on local machines, supported cloud providers, or Kubernetes clusters.

Key features include:
- Support for multiple backends (Docker, Kubernetes, AWS EC2, etc.)
- Based on the DevContainer standard
- Development environments as code
- Fast environment startup and teardown

## DevPod Provider and Workspace Creation

### Installing DevPod

DevPod offers a desktop version, but since I'm running it on FydeOS, I'll install the DevPod CLI. The DevPod CLI is quite friendly as it has almost no external dependencies and even includes built-in SSH functionality, making it easy to install on FydeOS's native command line.

```bash
curl -fsSL https://raw.githubusercontent.com/loft-sh/devpod/master/scripts/install.sh | sh
```

### Configuring Provider

Provider is the backend driver for DevPod, defining where the development environment runs. Using SSH as an example:

```bash
devpod provider add ssh --name amd -o HOST=AMD
```

The above addition method is based on existing SSH configuration, corresponding to the content in `~/.ssh/config`:

```
Host AMD
  HostName 192.168.0.2
  Port 22
  User user
  IdentityFile /home/user/.ssh/AMD_Key
```

Note: DevPod is a tool based on container technology, so the target machine for SSH needs to have Docker or Podman already installed.

Additionally, when using Ubuntu, pay special attention to the source of the Docker software. If Docker is installed via Snap on Ubuntu, due to Snap restrictions, it cannot access hidden folders (those starting with `.`) in the user's Home directory, causing DevPod to fail when executing any build commands with Docker.

### Creating Workspace

Specify the provider and project path to start the workspace. DevPod will create the necessary development environment based on the configuration files. More importantly, DevPod automatically installs the web version of OpenVSCode in the workspace, starts the service, automatically forwards the service port to the local machine, allowing you to code happily. This truly achieves the goal of development environments as code. At the same time, these configurations only require the server to have SSH port open, which is really quite convenient.

```bash
devpod up --provider amd --source git https://github.com/your-username/your-repo
```

This will create a development environment based on the configuration files in the `.devcontainer` folder of the repository.

## Introduction to DevContainer Format

DevContainer is a standard format for development environment configuration. This example includes two main files: `devcontainer.json` and `Dockerfile`.

### Dockerfile Configuration

The Dockerfile is essentially the configuration file used during docker build. DevPod uses this file to create the workspace container.

Note: It's best to build configuration containers based on containers from mcr.microsoft.com/devcontainers, as these containers are built according to specific requirements. Some of DevPod's features depend on files and configurations within the container. If using other containers, you may need to handle compatibility issues yourself when starting the workspace.

```dockerfile
# Use MS official base image, which contains vscode user to avoid issues when DevPod starts
FROM mcr.microsoft.com/devcontainers/base:bookworm

# Install pixi for vscode
USER vscode

# Install pixi, as my projects all use it for dependency management
RUN curl -fsSL https://pixi.sh/install.sh | sh

# Configure git user information and inject needed alias content into the environment
RUN git config --global user.name "Sylens" && \
  git config --global user.email "qiumin14@163.com" && \
  echo 'alias aider="aider --no-check-update --no-show-model-warnings --yes --no-auto-commits --model deepseek/deepseek-reasoner"' >> ~/.bashrc

ENV PATH="/home/vscode/.pixi/bin:${PATH}"
```

Another important point to note is that DevPod's operation logic is to first build the container, then mount the project code and perform further configuration according to `devcontainer.json`. Therefore, any settings that depend on project code cannot be performed during the container build phase.

### devcontainer.json Configuration

`devcontainer.json` defines the metadata and IDE configuration for the development container. Here's a basic configuration that only sets some environment variables and the minimum required plugins.

```json
{
  "name": "silen_blog",
  "build": { "dockerfile": "Dockerfile" },
  "postCreateCommand": "pixi run deploy",
  "containerEnv": {
    "AIDER_DARK_MODE": "${localEnv:AIDER_DARK_MODE}",
    "AIDER_CODE_THEME": "${localEnv:AIDER_CODE_THEME}",
    "DEEPSEEK_API_KEY": "${localEnv:DEEPSEEK_API_KEY}"
  },
  "customizations": {
    "vscode": {
      "settings": {
        "workbench.colorTheme": "Solarized Dark"
      },
      "extensions": [
        "naumovs.color-highlight",
        "ms-ceintl.vscode-language-pack-zh-hans",
        "tamasfe.even-better-toml"
      ]
    }
  }
}
```

Note: `postCreateCommand` can only contain one item. If you need to execute multiple commands, it's recommended to write them into a shell script, place it in the `.devcontainer` directory, and then call it with `postCreateCommand`.

## Afterword

Back when I was writing [ReviewGPT](https://github.com/SilenWang/ReviewGPT) last year, I realized that everyone has ideas, but what matters is actually building them.

I was originally modifying the terminal app in ChromeOS, hoping to add port forwarding and automatically install VSCode web on the target machine, then automatically forward the port to the local machine. As I worked on it, I discovered that Google had long abandoned pNaCl (which means native apps can't access local ports), so I thought about developing a Go application that meets my requirements. Continuing to work on it, I discovered DevPod - except that it's container-based... it's already very close to my needs... Well, what else should I write then... Let me first learn how to use DevPod...