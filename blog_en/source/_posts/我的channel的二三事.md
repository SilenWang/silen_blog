---
title: About My Channel
date: 2026-05-08 07:00:00
tags:
  - channel
  - pixi
  - GitHub Actions
categories: Others
---

I maintain a package channel on [prefix.dev](https://prefix.dev) to distribute some self-built packages. Since I've been busy with work recently, I haven't looked after it for a while. Today I checked and noticed that the Actions had mysteriously stopped running... So I looked into the reason and decided to document a few other issues as well.

<!-- more -->

## Quota Limits

Although prefix.dev's open channels are quite generous, they are not unlimited. Each channel has storage and package count quotas. Currently, my channel only has 5 packages but has already consumed 10GB of space. This is mainly because opencode updates very frequently — a single minor version can have dozens of sub-versions. Although each build is only a few dozen MB, the sheer number of versions adds up.

Currently, I have to manually check the channel for old packages and delete those that are no longer needed or have been superseded by newer versions. This process is entirely manual and quite tedious, but I haven't found a better automation solution yet.

## GitHub Actions Cron Jobs Are Automatically Disabled

I rely on GitHub Actions to automatically build and update packages in my channel. However, GitHub has a policy: if a repository's Actions have **no activity for 60 days**, the **scheduled workflows** are automatically disabled.

This means I need to manually re-enable these workflows on the GitHub repository page every two months or so. Without this, the channel's automatic updates would stop. So for now, the most worry-free approach is to set a reminder on my phone to re-enable the workflows every 60 days. It's not elegant, but it's simple and reliable.

## Channel Usage

At first I thought only I was using this channel, but it turns out others have been installing packages from it as well.

The most downloaded package is **opencode**. However, opencode is now also available via GitHub Releases for direct download, so fewer users should be installing through the channel going forward. This also means my channel maintenance pressure will gradually decrease — at least from the perspective of the opencode package.
