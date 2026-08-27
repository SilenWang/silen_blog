---
title: Netcatty - An SSH Client with Built-in AI
categories: Others
date: 2026-08-27 10:00:00
tags: ['Netcatty', 'SSH', 'SFTP', 'AI Agent', 'Open Source', 'Terminal']
---

I've been using multiple devices for a long time — my own laptop, servers, NAS, company servers. All of them accessed via SSH clients. Honestly, logging into servers manually and typing out commands for cleanup, maintenance, and troubleshooting gets pretty tedious.

To save myself the trouble, I've been using [opencode](https://opencode.ai/) to do cleanup and maintenance work directly on the target machines. It works, sure, but there's one unavoidable hassle: **you have to install opencode on that machine first**. You have to deal with environments, dependencies, keys, and if you also need to work inside containers, it gets even more annoying.

So I've been keeping an eye on tools like Wrap — SSH clients with built-in AI. But I didn't want to pay extra just for that. Eventually I found [Netcatty](https://github.com/binaricat/Netcatty).

<!-- more -->

## What is Netcatty

Netcatty is an **open-source SSH/SFTP client with built-in AI**, built with Electron, React, and xterm.js, supporting macOS, Windows, and Linux. On its own it's a fully-featured SSH workspace that covers most of the common needs:

- **Dual-pane SFTP browsing** with a built-in editor, so you can edit remote files right inside the client
- **Split terminals + tabs + session management**, so you can work on multiple servers side by side
- **Vault views** — organize hosts in grid, list, or tree layouts with fast search
- **Custom themes and keyword highlighting** — easier on the eyes, and clearer when debugging output

In other words, even ignoring its AI features, Netcatty is already a solid, modern SSH client.

## Built-in and custom agent support

Of course, the most important part is the built-in AI. Netcatty ships with an AI assistant called Catty Agent that can directly read the terminal output over SSH. Because of that, it understands your server environment, can execute commands, and even coordinate tasks across multiple hosts. For example, ask it to "check the load on this machine," and it runs the right commands, analyzes the output, and gives you a clear conclusion — instead of making you type `top` yourself and interpret it by hand.

Beyond the built-in Agent, it also supports **custom agents**. This matters to a lot of people, because most of us already have agents we're used to (Claude Code, Codex, etc.). Being able to plug in the setup you already know, rather than being locked into one specific agent, is clearly a plus.

## Core advantage

The advantage is obvious: **the agent doesn't need to be installed on the target machine.** The client establishes the SSH connection itself and reads the remote terminal through SSH, while the agent still runs locally. So all the environment, models, and workflows I've configured on my local machine work as-is, but they operate on servers thousands of miles away. That eliminates the hassle of setting up environments and dependencies over and over on each target machine.

For someone like me who uses multiple machines regularly, that difference is a real timesaver. The install-and-configure dance that used to repeat a dozen times now happens zero times.

## A few thoughts

A note off to the side: Netcatty certainly isn't the first SSH client with built-in AI (I don't even know if Wrap was the first), but among currently **fully open-source** options, it looks like the most polished one — features, cross-platform support, and AI integration are all quite complete. It also reinforces something I've been feeling more and more: many tools that have been around for a decade or more gain surprising convenience once they're combined with AI. Pain points that used to be "I need a missing command-line tool" now often resolve with a sentence in natural language plus AI.

If you also deal with a pile of remote servers and want to skip the "install the tool on every machine" hassle, Netcatty is worth a try. The repo is at [github.com/binaricat/Netcatty](https://github.com/binaricat/Netcatty), licensed under GPL-3.0. You can grab it, or check out [netcatty.app](https://netcatty.app) first.