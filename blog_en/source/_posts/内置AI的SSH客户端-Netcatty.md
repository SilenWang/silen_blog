---
title: Netcatty - An SSH Client with Built-in AI
categories: Others
date: 2026-08-27 10:00:00
tags: ['Netcatty', 'SSH', 'SFTP', 'AI Agent', 'Open Source', 'Terminal']
---

I maintain quite a few servers, and I've been using a traditional SSH client for the job. Honestly, logging into servers manually and typing out command after command for cleanup, maintenance, and troubleshooting gets old fast. Once AI became a thing, I started daydreaming about a tool where an agent could work directly on the servers while I just reviewed the results on the side.

I previously tried using [opencode](https://opencode.ai/) directly on the target machines for cleanup and maintenance work. The idea was sound, but there was one unavoidable hassle: **you have to install opencode on that machine first**. Setting up the environment, dependencies, and permissions for each one was such a chore that it ate away most of the convenience I was hoping for. So I always felt that approach was merely "workable," not "smooth."

It wasn't until I came across [Netcatty](https://github.com/binaricat/Netcatty) that I felt this kind of tool finally came close to what I wanted.

<!-- more -->

## What is Netcatty

Netcatty is an **open-source SSH/SFTP client with built-in AI**, built with Electron, React, and xterm.js, supporting macOS, Windows, and Linux. On its own it's a fully-featured SSH workspace that covers most of the common needs:

- **Dual-pane SFTP browsing** with a built-in editor, so you can edit remote files right inside the client
- **Split terminals + tabs + session management**, so you can work on multiple servers side by side
- **Vault views** — organize hosts in grid, list, or tree layouts with fast search
- **Custom themes and keyword highlighting** — easier on the eyes, and clearer when debugging output

In other words, even ignoring its AI features, Netcatty is already a solid, modern SSH client that works fine as a replacement for PuTTY, Termius, and the like.

## Built-in and custom agent support

What really caught my eye is its AI side. Netcatty ships with a built-in AI assistant called Catty Agent that understands your server environment, executes commands, and can even coordinate tasks across multiple hosts. For example, ask it to "check the load on this machine," and it runs the right commands, analyzes the output, and gives you a clear conclusion — instead of making you type `top` yourself and interpret it by hand.

But honestly, the more important thing for me is that it supports **custom agents**. That matters for my use case, because I already have an agent workflow I'm used to. Being able to plug in the setup I'm familiar with, rather than being locked into a built-in assistant, is a big plus in my book.

## What it actually solves

Thinking back to the pain of using opencode directly on servers, Netcatty's approach is quite different:

**The agent doesn't need to be installed on the target machine.** The client establishes the SSH connection itself and reads the remote terminal through SSH, while the agent still runs locally. So all the environment, models, and workflows I've configured on my local machine work as-is, but they operate on servers thousands of miles away. That eliminates the hassle of setting up environments and dependencies over and over on each target machine.

For someone like me who maintains a fleet of machines, that difference is a real timesaver. The install-and-configure dance that used to repeat a dozen times now happens zero times.

## A few thoughts

A note off to the side: Netcatty certainly isn't the first SSH client with built-in AI, but among currently **fully open-source** options, it looks like the most polished one — features, cross-platform support, and AI integration are all quite complete. It also reinforces something I've been feeling more and more: many tools that have been around for a decade or more gain surprising convenience once they're combined with AI. Pain points that used to be "I need a missing command-line tool" now often resolve with a single sentence in natural language.

If you also deal with a pile of remote servers and want to skip the "install the tool on every machine" hassle, Netcatty is worth a try. The repo is at [github.com/binaricat/Netcatty](https://github.com/binaricat/Netcatty), licensed under GPL-3.0. You can grab it, or check out [netcatty.app](https://netcatty.app) first.