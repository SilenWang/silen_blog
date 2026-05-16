---
title: Tales of My Channel
date: 2026-05-08 07:00:00
tags:
  - channel
  - conda
  - GitHub Actions
categories:
  - Tech Talk
---

## Preface

I maintain a conda channel on [prefix.dev](https://prefix.dev) to distribute some self-built packages. Over time, the maintenance workload has gradually exceeded my expectations. Here I'd like to document some of the issues I've encountered.

<!-- more -->

## Quota Limits

Resources on prefix.dev channels are not unlimited. Each channel has storage and package count quotas. As I uploaded more and more packages to my channel, I quickly hit the quota ceiling.

Currently, I periodically check the channel for old packages and manually delete those that are no longer needed or have been superseded by newer versions. This process is entirely manual and quite tedious, but I haven't found a better automation solution yet.

## The GitHub Actions Cron Job Dilemma

I rely on GitHub Actions to automatically build and update packages in my channel. However, GitHub has a policy: if a repository's Actions have **no activity for 60 days**, the scheduled workflows are automatically disabled.

This means I need to manually re-enable these workflows on the GitHub repository page every two months or so. Without this, the channel's automatic updates would stop.

I've tried a few workarounds:
- Adding a dummy commit to trigger the workflow? — This is still a manual operation at its core.
- Using an external service to periodically ping the repo? — This adds extra dependencies and complexity.

For now, the simplest approach is to set a reminder on my phone to re-enable the workflows every 60 days. It's not elegant, but it's simple and reliable.

## Channel Usage

At first I thought only I was using this channel, but it turns out others have been installing packages from it as well.

The most downloaded package is **opencode**. However, opencode is now also available via GitHub Releases for direct download, so fewer users should be installing through the channel going forward. This also means my channel maintenance pressure will gradually decrease — at least from the perspective of the opencode package.

## Afterword

Maintaining a conda channel comes with its share of tedious hassles, but overall it's still been a valuable experience. If you're considering setting up your own channel, I hope this little article helps you anticipate some of the pitfalls you might encounter.

If you also maintain a conda channel, feel free to share your experiences!
