---
title: Monitoring My Websites with Uptime Kuma
categories: Coding
date: 2025-10-09 15:30:29
tags: [Monitoring, Uptime-Kuma, Website Maintenance, Service Monitoring]
---

Uptime Kuma is an extremely user-friendly open-source monitoring tool that makes monitoring various network services simple and efficient. It supports multiple protocol monitoring including HTTP(s), TCP, Ping, DNS queries, and even simulates real user access through Chrome engine mode to more accurately monitor service stability. Additionally, it supports dozens of notification methods for alerts and can integrate with various applications we use daily. It also offers a unique Push monitoring method, allowing users to extend monitoring functionality through custom code (such as monitoring SSH service availability).

<!-- more -->

## Introduction to Uptime Kuma

[Uptime Kuma](https://github.com/louislam/uptime-kuma) is a service monitoring application developed by Louis Lam. The design philosophy of this project is simplicity and ease of use, making it accessible even for users without technical backgrounds.

## Installing Uptime Kuma

The simplest method is to use Docker for installation. In iStoreOS's app store, there are pre-configured solutions available - just install and start it. The project runs on port `3001` by default. If you're also using Gogs, there might be a port conflict, so you'll need to change the port.

## Basic Usage

### 1. Initial Setup

When accessing for the first time, you need to set up an administrator account and password. The entire process takes only a few minutes to complete.

### 2. Adding Monitoring Items

Click the "Add Monitor" button and fill in the following information:

- **Monitor Type**: Select the monitoring type (such as HTTP, TCP, Ping, etc.)
- **Friendly Name**: Give the monitor an easy-to-remember name
- **URL/Hostname**: The website address or host to monitor
- **Heartbeat Interval**: Check interval (default 60 seconds)

Particularly worth mentioning is the **Chrome Engine Monitoring Mode**: In the advanced settings of HTTP monitoring. This mode uses the actual Chrome browser kernel to access websites, can execute JavaScript, and detect whether page loading is truly complete. This is especially useful for Single Page Applications (SPA) and websites that require user interaction.

### 3. Setting Up Notifications

In the "Setup Notification" section, you can add various notification methods. Uptime Kuma supports dozens of notification channels, including:

- **Instant Messaging**: Telegram, Discord, Slack, Line, Mattermost, etc.
- **Email**: SMTP email notifications
- **Mobile Push**: Pushover, Gotify, NTFY, etc.
- **Webhook**: Can be configured to directly integrate with domestic office software like WeCom, DingTalk, Feishu, etc.

Taking WeCom as an example:
1. Create a notification bot in a group chat
2. Note down the bot's key
3. Select WeCom notification method in Uptime Kuma
4. Fill in the key

This way, when services become unavailable or recover, you'll receive notifications in WeCom.

### 4. Using Push Monitoring Mode

Push monitoring is a unique feature of Uptime Kuma that allows services to actively report their status to Uptime Kuma. This means even for unsupported services (like SSH), we can write simple scripts to push monitoring results to Uptime Kuma for monitoring. I plan to write a Go service to monitor SSH connectivity when I have time.

## Summary

Uptime Kuma is an extremely feature-rich and very user-friendly monitoring tool, particularly suitable for individual developers and small teams.

If you're also looking for a powerful yet easy-to-use monitoring solution, I highly recommend trying Uptime Kuma - let's become cyber monitoring room masters together!
