---
title: Compile openfyde
categories: Others
date: 2025-01-04 14:58:33
tags: ['openfyde', 'fydeos', 'chromeos']
---


Compile the OpenFyde Image for Fydetab Duo by synthesizing official instructions and forum articles, aiming to compile the latest R120 version image.

<!-- 摘要部分 -->
<!-- more -->

Fyde OS originally provided a domestic accelerated mirror, but it seems to be unavailable now. Therefore, a VPN is required to complete the following content.

## System Preparation

- System: Ubuntu 22.04, AMD Ryzen 5 1600 Six-Core @ 12x 3.2GHz, 64G
- Install necessary software: `sudo apt-get install git gitk git-gui curl xz-utils python3-virtualenv python3-oauth2client lvm2 thin-provisioning-tools repo`
- Install `depot_tools`:
    + Clone the project: `git clone https://chromium.googlesource.com/chromium/tools/depot_tools.git && cd depot_tools && git checkout e121d14b12412e95ac833cfd31602b674499ea25`
    + Set environment variables: `export PATH=/mnt/hdd1/chromeos/openfyde/depot_tools:$PATH`
    + You can also set `export DEPOT_TOOLS_UPDATE=0` to avoid tool updates. If the network is stable, this setting is optional.

## Api Key Preparation
- Since I am still using a Google account under Fydeos, I applied for a Google API key instead of Fydeos's. Follow the steps in the [official documentation](https://www.chromium.org/developers/how-tos/api-keys/)
- General steps:
    + Join the chromium-dev development group
    + Login to `https://cloud.google.com/console`
    + Create a new project, e.g., openfyde.
    + Select API Services > Library, find the following API projects, and enable them one by one:
        * Cloud Search API
        * Google Drive API
        * Safe Browsing API
        * Time Zone API
    + Go to API Services > Credentials, configure OAuth information as prompted, then create credentials > Create OAuth client ID. I chose the desktop application type (no Others option as described in the documentation).
    + Create credentials > Create API key to obtain the key.
    + Save the client ID, client secret, and API key, and store them in the `~/.googleapikeys` file as follows:
        ```txt
        'google_api_key': 'your api key',
        'google_default_client_id': 'your client id',
        'google_default_client_secret': 'your client secret',
        ```

## Source Code Clone

- Create a directory: `mkdir -p r120/openfyde`
- Configure Git, otherwise cloning the code will result in errors:
    ```bash
    git config --global user.name "sylens"
    git config --global user.email "silenseek14@gmail.com"
    ```
- Clone the source code: The R126 version of `openFyde/manifest.git` is not released，so only R120 can be compiled. `repo sync -j8` requires downloading a large amount of source code, with chromium being the largest at nearly 40G. If the download fails due to network issues, running it again will restart the download. Ensure your VPN has sufficient bandwidth.
    ```bash
    repo init -u https://chromium.googlesource.com/chromiumos/manifest.git --repo-url https://chromium.googlesource.com/external/repo.git -b release-R120-15662.B

    git clone https://github.com/openFyde/manifest.git openfyde/manifest -b r120-dev
    ln -snfr openfyde/manifest .repo/local_manifests
    repo sync -j8
    cd openfyde/chromium
    gclient sync
    ```

## Entering the Chroot Environment

- Return to the `r120` directory and run `cros_sdk` to start downloading the SDK and preparing the environment. This process takes a long time.

## Starting the Compilation

- Since the current version of archero is unusable, modifications are needed to avoid errors.
- After entering chroot, you will be in `/mnt/host/source/src/scripts` by default. The files that need modification are:
    + `/mnt/host/source/src/overlays/overlay-fydetab_duo-openfyde/metadata/layout.conf`: Delete `archero`、`tablet`、`ai-dev` in the first line
    + `/mnt/host/source/openfyde/overlays/overlay-fydetab_duo-openfyde/profiles/base/parent`: Delete the second line's`archero:base`
- Note: All compilation commands are executed within chroot.
- Run `setup_board --board=fydetab_duo-openfyde`
    + It appears that the setup is not fully completed, but the target directory has been generated.
- Install `capnproto`: `sudo emerge capnproto`
- Modify `/mnt/host/source/chromite/lib/dlc_allowlist.py`: Change line 14 to `DLC_FACTORY_INSTALL = (r"termina-dlc", r"sample-dlc",)`, i.e., add `r"termina-dlc"`
- Run `cros build-packages --jobs=4 '--board=fydetab_duo-openfyde' --no-withautotest --autosetgov --no-use-any-chrome`
- After completion, build the image with `cros build-image --board=fydetab_duo-openfyde --noenable_rootfs_verification`

## Other Issues During Compilation

- Compilation interrupted abnormally: Encountered an abnormal interruption once. Re-entering the environment and compiling again caused inexplicable issues. To start over, delete the `chroot` and `out` folders, then run `cros_sdk` again from the beginning.

- Virtual hard disk settings: As [hinted by a forum expert](https://community.fydeos.com/t/topic/48326/14), an additional setting is required:
    + Open `/etc/sysctl.conf` and add `vm.max_map_count=524288`.
    + Run `sudo sysctl -p`.
