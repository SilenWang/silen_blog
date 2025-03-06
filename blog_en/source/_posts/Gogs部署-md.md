---
title: Gogs Deployment
date: 2018-09-21 18:32:00
tags: [gogs, vps]
categories: Others
---

Bought a VPS; the monthly traffic usage is so low that it's really a waste... So let's make better use of it. After all, it costs $5 per month.

<!-- more -->

# Deployment Steps
Gogs is a lightweight Git hosting service implemented in Go. Its resource consumption is so small that it can even run on a Raspberry Pi.

I previously set it up once but didn't document the process. This time, I had to look it up again...

Overall, using Docker for deployment (it's simple), so basically just configuring a file and running a few commands. The installation of Docker isn't explained here; there are instructions in ArchWiki. The main installation command is:

```bash
docker pull docker.io/gogs/gogs
```

Then prepare a `docker-compose.yml` file with the following content, execute it in the directory containing this file, and start using Gogs. The access port for Gogs is 10080.

```txt
version: "2"

services:
  gogs:
    image: docker.io/gogs/gogs

    restart: always
    ports:
      - "10022:22"
      - "10080:3000"
    volumes:
      - /var/gogs:/data

```

Then to enable passwordless push, you need to set up SSH. By default, Gogs doesn't have the SSH server enabled. You need to edit `/var/gogs/gogs/conf/app.ini` and modify it as follows (this is my configuration file directory; adjust according to the actual mount path in `volumes`):

```ini
[server]
DOMAIN           = http://XXX.XXX.XXX.XXX
HTTP_PORT        = 3000
ROOT_URL         = http://localhost:10080/
DISABLE_SSH      = false
SSH_PORT         = 10022
START_SSH_SERVER = true
OFFLINE_MODE     = false
```

Then set up the SSH keys as you would with GitHub.

# Solving the Issue of Adding Multiple Keys (update@20180921)

Although Gogs has a graphical interface for setting up SSH keys, there might be a bug in the version I pulled. When adding the first key, it works fine, but when trying to add more than one key from multiple machines, even though it shows as added in the graphical interface, it can't normally log in via SSH with passwordless access. You need to further set up:

- Since my Docker container mounts `/data` to a directory outside the container, there's no need to log into the Docker system and modify directly.
- The corresponding file is: `挂载目录/gogs/git/.ssh/authorized_keys`. Add public key content according to the format. For example:
```txt
command="/app/gogs/gogs serv key-1 --config='/data/gogs/conf/app.ini'",no-port-forwarding,no-X11-forwarding,no-agent-forwarding,no-pty ssh-rsa AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA test
```

- Save it and you can log in without restarting the container.
