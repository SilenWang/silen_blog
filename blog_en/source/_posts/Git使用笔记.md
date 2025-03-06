---
title: Git Usage Notes
categories: Others
date: 2018-09-28 12:34:12
tags: ['git']
---

I had read about Git tutorials a long time ago, but I never actually used it. Now that I'm using it, I feel like I've become a programmer in some way!

# Operation Method Records

- User Initialization, when first using on a new device:

```bash
git config --global user.name "Silen Wang" # Set the user name
git config --global user.email "mymail@gmail.com" # Set the user email
git config --global core.editor 'vim' # Set the editor used for commit messages
```

- Clone a project's entire code:

```bash
git clone PROJ_URL
```

- Pull updates

```bash
git pull # Applicable when cloning a project, if there are multiple sources, specify which source to pull from
git pull origin master # Specify pulling the master branch of the origin remote repository
```

- Push updates

```bash
git push
git push origin master # Specify pushing to the master branch of origin, if not specified, it might use a default setting? Not sure exactly what happens.
```

- Remote management

```bash
git remote # View remote information
git remote add REMOTE_URL # Add remote repository information
```

- Modify commits, Git stores local file modifications in a staging area. To submit these changes, first confirm the commit modifications.

```bash
git add FILE # Stage new files or mark changes to files as ready for submission
git commit # Commit all uncommitted modifications with a simple description of what was modified. This description can be used to trace back previous actions.
git commit -a # Directly commit all uncommitted modifications with the same description
```

- Branch operations, the concept of branches allows creating a backup of the current project (a copy), then making modifications on this branch without affecting the original branch (e.g., master). Once the modified branch is tested and has no issues, merge the changes into the branch that needs modification.

```bash
git checkout -b dev # Create and switch to the dev branch
git checkout dev # Switch to the dev branch
git merge dev # Merge changes from the dev branch into the current branch
```
```