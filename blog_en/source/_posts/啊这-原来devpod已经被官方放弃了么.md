---
title: Wait... Has Devpod Really Been Abandoned by Its Maintainers?
categories: Others
date: 2026-01-17 02:54:00
tags: 
  - Devpod
  - Open Source Maintenance
  - Community Fork
  - Project Sustainability
  - Docker Compose
  - Open Source Software
---

It all started when I tried to use Devpod to set up a container based on Docker Compose, only to find that Devpod didn't seem to be able to call Docker Compose correctly to create the container. I didn't think much of it until I searched online... Wait, what? Has this project really been dropped by loft‑sh?

<!-- more -->

According to [issue 1946](https://github.com/loft-sh/devpod/issues/1946) and [issue 1915](https://github.com/loft-sh/devpod/issues/1915), vCluster is loft‑sh’s most commercially successful and promising project, so all their energy is temporarily focused on it, leaving no time to maintain community projects.

This situation feels similar to the Fyde team—they also mentioned on Discord that they would prioritize work that keeps the company running...

Fortunately, a [community hero](https://github.com/skevetter) has created a community fork of Devpod and is updating it at an impressive pace (already at version 0.9.*)!

Using skevetter’s fork, Docker Compose works normally, which is a silver lining in an otherwise unfortunate situation.

However, so far, the vast majority of the code has been contributed by skevetter alone, which isn’t really healthy. Just like with onelist before, the project became quite popular, but almost all development was concentrated on the original author. Relying on pure passion isn’t sustainable in the long run and partly led to the eventual sale and abandonment of that project.

I hope Devpod can follow the path of Gitea, which originated from Gogs but eventually grew a complete community, allowing a good project to continue thriving.
