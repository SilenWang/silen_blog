---
title: Deploying RStudio with Singularity and a Conda Environment
tags: ['conda', 'rstudio', 'singularity']
categories: Bioinfomatic
date: 2023-07-26 16:13:50
---

Today, I needed to deploy an RStudio server for someone else. Here are the steps and key points of the deployment.
<!-- Abstract part -->
<!-- more -->

## Introduction
Running RStudio requires root permissions and cannot be installed directly using conda or mamba like Jupyter. Therefore, this deployment used Singularity, which is a container technology I had never used before. While Sigularity differs from Docker or Podman in design, it is primarily designed for HPC environments with different permission management. However, these differences do not significantly affect its use in this case.

## 1. Container Preparation
Although Sigularity and Docker are different container technologies, Sigularity provides comprehensive compatibility support for Docker images. Built Docker containers can be easily converted to the `sif` format used by Singularity. For this deployment, I referenced the instructions from the [rstudio-server-conda](https://github.com/grst/rstudio-server-conda/blob/master/README.md) project and pulled the Docker image using `singularity pull docker://rocker/rstudio:4.2`.

## 2. Necessary File Preparation

I used the following bash script to generate several directories and configuration files.

```bash
mkdir -p run var-lib-rstudio-server
printf 'provider=sqlite\ndirectory=/var/lib/rstudio-server\n' > database.conf
echo "auth-minimum-user-id=100" > rserver.conf
echo "session-default-working-dir=/home/silen/Rstudio/Workspace" > rsession.conf
echo "session-default-new-project-dir=/home/silen/Rstudio/Workspace" >> rsession.conf
```

This script generates `run` and `var-lib-rstudio-server` directories, as well as `database.conf`, `rserver.conf`, and `rsession.conf` configuration files. The `auth-minimum-user-id=100` setting in `rserver.conf` is crucial because the RStudio service runs under the rstudio user with UID 999 in the container. If this user does not exist or if the permissions are not switched to this user, the container will fail to run. I referenced a post on GitHub ([#18](https://github.com/grst/rstudio-server-conda/pull/18)) and used the `--server-user` parameter to specify that the service should be run using the current user.

However, by coincidence, my current user's UID is less than 513, while the default `auth-minimum-user-id` value is 1000, which prevents users with UIDs less than 1000 from running the service (apparently for security reasons). Therefore, I needed to change the parameter value in `rserver.conf` to allow the service to start normally.

## 3. Starting the Container with Specified Parameters

I used the following command to start the service. The main commands bind the previously created directories and configuration files to the container so that the service can read them. Additionally, I specified the IP address and port for easy access. Since this is a shared use scenario, I did not set up authentication. I also referenced the [rstudio-server-conda](https://github.com/grst/rstudio-server-conda/blob/master/README.md) project to specify the conda environment used by RStudio, so there's no need to reinstall R packages within the container.

```bash
singularity exec \
   --bind run:/run \
   --bind var-lib-rstudio-server:/var/lib/rstudio-server \
   --bind database.conf:/etc/rstudio/database.conf \
   --bind rserver.conf:/etc/rstudio/rserver.conf \
   --bind rsession.conf:/etc/rstudio/rsession.conf \
   --bind /home/silen/R-Plot:/etc/R-Plot \
   --env CONDA_PREFIX=/etc/R-Plot  \
        --env RSTUDIO_WHICH_R=/etc/R-Plot/bin/R \
   rstudio_4.2.sif \
   /usr/lib/rstudio-server/bin/rserver --www-address=0.0.0.0 --www-port=7788 --server-user=$(whoami)
```

After starting, the interface looks like this:

![ui](https://raw.githubusercontent.com/SilenWang/Gallary/master/2023/07/upgit_20230727_1690390228.png)

## 4. Notes

Sigularity feels more similar to snap and flathub software packaging projects, where the host and container are not completely isolated from each other. After starting RStudio, it can directly access the user's Home directory on the host machine. I'm not sure if this is a special behavior of this container or if Sigularity behaves similarly in general. We will need to do more dependency analysis and containerization for our internal analysis processes... It seems we still need to learn more about the principles and usage of Singularity.
```