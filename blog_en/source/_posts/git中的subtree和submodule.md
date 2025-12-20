---
title: Subtree vs Submodule in Git
categories: Other
date: 2025-12-06 01:49:33
tags:
  - git
  - version control
  - subtree
  - submodule
---

Emmmmm, the contract with our official website vendor expired this year, so... I ended up with one more project to manage... This project is still frontend/backend separated, but the difference is that the official site has an English version, and the English version is actually a branch of the frontend project. I’ve used Submodule before to integrate a colleague’s independent module into the main project, but this time, after consulting with AI, I chose a different approach: Subtree.

<!-- more -->

## Quick Comparison: Submodule vs Subtree

In short, Submodule creates a “link” in the main repository that points to a specific commit of a sub‑repository, while Subtree copies the entire code of the sub‑repository into a sub‑directory of the main repository. Both have their pros and cons; the table below helps you quickly grasp their core differences:

| Feature | Submodule | Subtree |
|---------|-----------|---------|
| How code is stored | Main repo only keeps a reference (commit hash) to the sub‑repo | Sub‑repo’s code is fully copied into the main repo |
| Extra initialization after cloning | Requires `git submodule init && git submodule update` | No extra steps, all files are ready |
| Directly modifying sub‑project code inside the main project | Must enter the sub‑module directory and commit separately, which is cumbersome | Can modify and commit directly within the main repo directory |
| History | Sub‑repo history is completely separate from the main repo | Sub‑repo history can be optionally merged (`--squash`) or kept intact |
| Suitable scenarios | Sub‑project is relatively stable, updated infrequently; multiple main projects share the same sub‑project | Need to frequently modify the sub‑project inside the main project; need to bring in multiple different branches of the same repository |

In a nutshell, `Subtree` essentially clones a full copy of the sub‑project code. When you make changes and commit, the modifications are synchronized back to the sub‑project, and when the sub‑project is updated elsewhere, those changes can be pulled back into the main project. On the other hand, `submodule` is more like a soft link—the main project only specifies which version of the sub‑project code to use, and pushing changes back to the sub‑project is relatively more cumbersome.

For my case, on one hand I wanted **to bring two different branches (`master` and `en`) of the same repository as separate sub‑directories**; on the other hand, I actually need to modify the sub‑project code frequently, and the main project is essentially an integration of the frontend and backend. According to AI’s advice, **Subtree** is the simpler combination to operate.

## Practical Usage of Subtree

Below are the specific commands I used to add each sub‑project to the main project. Note that I used the `--squash` parameter, which merges the sub‑project’s historical commits into a single commit, preventing the main project’s history from being flooded with numerous sub‑project commits.

### 1. Add remote repositories
```bash
git remote add Backend https://github.com/user/Backend
git remote add Frontend https://github.com/user/Frontend
```

### 2. Bring sub‑projects in as subtrees
```bash
# Chinese site (master branch)
git subtree add --prefix=app/Backend Backend master --squash
# English site independent branch (en branch)
git subtree add --prefix=app/Backend_EN Backend en --squash

# Backend admin project
git subtree add --prefix=app/Backend Backend_Admin aiming_med --squash
```

### 3. Daily sync operations
- **Pull updates from the sub‑project**: `git subtree pull --prefix=<directory> <remote> <branch> --squash`
- **Push changes made to the sub‑project inside the main project**: `git subtree push --prefix=<directory> <remote> <branch>`

Of course, to simplify operations, you can configure corresponding commands in Pixi.

### 4. Notes
- After setting up subtree as described above, the sub‑project code already has a full copy in this project. Therefore, when deploying again, after cloning the main project you **do not** need to run `subtree add` again.
- If you clone for development purposes, you still need to add the remote repositories (as in step 1 above) so that you can later push updates to the sub‑project.

## Afterword

I’ve been using Subtree for two weeks now. Apart from stumbling over the fact that the git in Devcontainer didn’t have the subtree module, everything else has been working fine... Hopefully no new issues will pop up...
