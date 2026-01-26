---
title: Cherry-pick is Everywhere
categories: Coding
date: 2026-01-20 00:05:43
tags:
  - Git
  - Version Control
  - "cherry-pick"
---

I wanted to contribute to the Fydetab Duo Wiki, but to preview the blog locally, I needed to add pixi or other configurations to the project. These changes were not suitable for the original repository. So I learned how to selectively submit commits, which is where `git cherry-pick` comes in.

<!-- more -->

## Git's cherry-pick

In the Git version control system, the `cherry-pick` command is used to apply a specific commit to the current branch without merging the entire branch. This perfectly fits my needs: I forked the source code, created a branch (dev) for modifications, where the first two commits were files for local configuration, and the later ones were actual changes to the documentation. When I finished the modifications and was ready to submit, I could create a dedicated submission branch based on the unmodified original branch, then use `cherry-pick` to select the last three commits to include:

```bash
# Apply specified commits to the current branch
git cherry-pick COMMIT3 COMMIT4 COMMIT5

# If conflicts occur, resolve them and continue
git cherry-pick --continue

# Or abort the cherry-pick
git cherry-pick --abort
```

Then, using this branch to submit a PR ensures that unnecessary content isn't pushed to the original project.

At the same time, after the upstream merges my PR, I can pull the changes from upstream into the dev branch, allowing me to continue updating, and then modify and use the dedicated merge branch again.

## Cherry-picking in data analysis

When I first saw `cherry-pick`, I found it quite interesting because I had also encountered this term in data‑analysis discussions, where it carries a more negative connotation—referring to "selective presentation." Honestly, I've done quite a bit of this after starting work... It's one of those "charming" moments in statistics...

Also, I asked an AI, and it seems the phrase "cherry pick" has been around for a long time... Hmmm... Indeed, there's nothing new under the sun...
