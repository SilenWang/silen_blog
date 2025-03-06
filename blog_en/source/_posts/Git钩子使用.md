---
title: Git Hooks Usage
categories: Script
date: 2019-08-17 19:33:27
tags: ['git', 'githook']
---

Although I also use Git for some project code management, I actually don't know many of Git's features and haven't used them. Today, I tried using Git hooks.
<!-- Abstract part -->
<!-- more -->
Hook, directly translated as a hook, refers to a series of predefined scripts. These scripts are automatically executed when we use specific Git commands to perform some code management operations, with the purpose of reducing repetitive operations and improving efficiency.

For example, my blog source code is managed by Git. The steps I take to write a blog post and publish it are:
1. `hexo new post` creates a markdown draft for a new blog post.
2. Open the newly created markdown file and write content.
3. After writing, run `hexo g` and `hexo d` to generate and push the blog content to the Gitpages source.
4. Run `git commit` to submit changes to the blog source code, then `git push` to the source repository.

Among these steps, there may be some variable operations during the writing process, but the last two steps are always repeated and not very changeable each time I write a blog post. These steps can be executed using hooks.

First, we need to determine when these things should be executed: My blog source code is checked before `commit`. Therefore, theoretically, any hook triggered by `git commit` can be used. The actual selected hook is `post-commit`, which executes after `commit`, so I wrote the following content into the `.git/hooks/post-commit` file of the project (the name should not be changed; specific names are required to be executed).

```bash
#!/bin/bash
cd /home/silen/git_proj/silen_blog/
hexo g
hexo d
```

Then, give this file executable permissions. This way, every time I complete a commit for this project, the script will automatically generate and push the blog.

Of course, the example here may not be very attractive because these two steps are actually not many... saving a lot of effort. But now, in addition to my main site, I have set up a personal quick reference manual on a sub-site using mkdocs, which is another independent project. If I want to achieve that one update synchronizes the other, hooks would be very useful.

Finally, it's worth noting that Git hooks are divided into client-side and server-side. I am currently using a client-side hook. One drawback of client-side hooks is that they cannot synchronize across different copies of the project (different copies are different clients), which is a problem I will need to solve later...

That's all.
