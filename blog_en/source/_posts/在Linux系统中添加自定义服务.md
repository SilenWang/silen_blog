---
title: Adding Custom Services in Linux Distributions
categories: Others
date: 2019-10-05 16:00:31
tags: ['systemd']
---

When using a Linux system, it's often necessary to set up things that start automatically at boot. For distributions that use systemd, writing your own service files and enabling them is a great choice.

<!-- 摘要部分 -->
<!-- more -->

Service files are written in a configuration file format similar to ini. A basic service configuration file includes the following content:

- `Unit`: The defined section, such as what the service does
- `Service`: The actual execution content, including what type of service it is, what to do when starting, what to do when stopping, and what to do when restarting
- `Install`: It's mostly about dependencies, such as when the service should start or before or after XX starts

Here is an example of a v2ray service:

```ini
[Unit]
Description=Daemon to start V2ray

[Service]
Type=simple
ExecStart=/usr/bin/v2ray -config /home/silen/script/v2ray/v2ray.json

[Install]
Alias=v2ray
WantedBy=default.target
```

If this service is enabled, it will start automatically by default(`WantedBy=default.target`), executing the command`/usr/bin/v2ray -config /home/silen/script/v2ray/v2ray.json`, when booting, which means activating the v2ray proxy. Other actions for the service, such as restarting and stopping, are not defined.

After writing the above file, you can name the file`v2ray.service`and place it in the specified location (`/etc/systemd/system` or `/usr/lib/systemd/system`), then use `systemctl enable v2ray.service`with root privileges to activate the service (it will start according to the configuration after boot), and use `systemctl start v2ray.service` to start the service immediately.

It should be noted that when calling root privileges to start the service, the executing user is naturally root, which is not always a good choice. Therefore, my previous synchronization service and this v2ray service are actually both started with personal user permissions. In other words, put the above service file in `/home/<user>/.config/systemd/user`, and then directly `systemctl start --uesr v2ray.service` to start the service.

In this case, the service will only start after my current user logs in, and there is no need for root privileges (whether the user needs to be in the root group or not, I am not sure).

The v2ray program used in the example belongs to the type that occupies the foreground, so it can be directly used with the `sample` and `oneshot` type. In addition to this type, there are: `simple`, `forking`, `oneshot`, `notify`, `dbus`, `idle`.

Of course, my custom services are generally quite simple, `simple` and `oneshot` should be enough.

As for `ExecStart`, `ExecStop`, and `ExecReload`, the former is required, while the latter two can be written or not depending on the start command situation.

Finally, `WantedBy=default.target` in `Install` is a must according to my test, because the service will not start automatically without it... It's not clear whether this is only for services at the user level.