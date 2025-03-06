---
title: Using Syncthing for File Synchronization
categories: Others
date: 2019-08-18 00:03:33
tags: ['Syncthing']
---

As long as you use multiple devices for work, you will inevitably encounter the problem of mutual synchronization between different devices. Now that it's not the era when cloud services just started to rise and were cheap and generous, there probably isn't a reliable and user-friendly service provider that won't suddenly run away or change its terms of service. So, we still have to rely on ourselves to set up our own...
<!-- Abstract part -->
<!-- more -->

Syncthing is introduced in [Little Known Software](https://www.appinn.com/syncthing/) and [Anycard Software](https://www.iplaysoft.com/syncthing.html), a simple yet open-source synchronization software.

Since the official version already provides pre-packaged executable files, Syncthing can be run directly on three major platforms without needing to download any additional packages. Therefore, I won't repeat that part here. This post mainly records how to set up Syncthing as a service on my computer and VPS.

Both CentOS7 and Manjaro use systemd, and the command for managing services is `systemctl`, so the operations are the same.

First, find the Syncthing executable file package. Inside the `Syncthing/etc/` directory, there are sample service configuration files. I used the `Syncthing/etc/linux-systemd/user/syncthing.service` file. Open this file and change the program path in the `[Service]` section of `ExecStart`. Since I am setting up this service as a user-level service rather than a system-wide service, I copied this file to `~/.config/systemd/user/`, then started the service using the following commands:

```bash
systemctl --user enable syncthing.service
systemctl --user start syncthing.service
```

After starting, you can use `systemctl` like other services to check the status of the Syncthing service.

![Syncthing](https://raw.githubusercontent.com/SilenWang/Gallary/master/Syncthing.png)
```