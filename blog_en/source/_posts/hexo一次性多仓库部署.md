---
title: Hexo One-Time Multi-Repository Deployment
categories: Others
date: 2020-08-02 01:19:38
tags: ['hexo', 'gitee']
---

Recently, GitHub's access has been increasingly problematic. Sometimes it's even difficult to check my blog posts... So I also have a copy hosted on the domestic Git hosting service Gitee. However, this raises a question: how can I update multiple places at once?

Firstly, set up the repository on Gitee. The basic method is similar to GitHub. The only thing you need to note is that Gitee's default project display is different from GitHub. GitHub uses `USERNAME.github.io`, while Gitee uses `USERNAME`.

Since Gitee has a feature to fork GitHub projects, when creating the project, choose to fork the original GitHub blog repository and then change the repository name to `USERNAME` (as shown in the example after migration... so it displays that the project has been created).

![gitee_proj](https://raw.githubusercontent.com/SilenWang/Gallary/master/gitee_proj.png)

Then, on the project page, select the service menu to set up access from `USERNAME.gitee.io` to view the blog.

![gitee_pages](https://raw.githubusercontent.com/SilenWang/Gallary/master/gitee_pages.jpg)

Next, set up SSH access to Gitee. This is similar to GitHub and won't be expanded here.

Then, configure the `_config.yaml` file in the project. In the deployment section, add the following configuration so that you can deploy twice at once.

```yaml
deploy:
  - type: git
    repo:
      git@github.com:SilenWang/silenwang.github.io.git
    branch: master
  - type: git
    repo:
      git@gitee.com:silenwang/silenwang.git
    branch: master
```

Originally, this should have been it... but when testing, another issue arose. My home network can no longer directly connect to GitHub, so the deployment to GitHub fails...

Therefore, I searched for ways to add a proxy to Git.

Since my configuration uses SSH protocol to connect to GitHub, after checking, it's necessary to add a `ProxyCommand` setting in the SSH configuration file as follows:

```config
Host github.com
    User git
    HostName github.com
    ProxyCommand nc -x 127.0.0.1:1080 %h %p
```

After adding this, I tried `hexo d` again and it successfully deployed.

That's all for now~
