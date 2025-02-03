---
title: 编译openfyde
categories: Others
date: 2025-01-04 14:58:33
tags: ['openfyde', 'fydeos', 'chromeos']
---


综合官方说明、论坛文章的内容，整理编译Fydetab duo的openfyde镜像的方法, 目标是编译最新的R120版本镜像
<!-- 摘要部分 -->
<!-- more -->

Fyde OS官方原本准备了国内的加速镜像，但是现在似乎已经不可用了，因此完成以下内容需要自备梯子

## 系统准备

- 编译系统为Ubuntu 22.04, AMD Ryzen 5 1600 Six-Core @ 12x 3.2GHz, 64G
- 安装必要软件: `sudo apt-get install git gitk git-gui curl xz-utils python3-virtualenv python3-oauth2client lvm2 thin-provisioning-tools repo`
- 安装`depot_tools`:
    + 克隆项目: `git clone https://chromium.googlesource.com/chromium/tools/depot_tools.git && cd depot_tools && git checkout e121d14b12412e95ac833cfd31602b674499ea25`
    + 设置环境变量: `export PATH=/mnt/hdd1/chromeos/openfyde/depot_tools:$PATH`
    + 可以同时设置`export DEPOT_TOOLS_UPDATE=0`避免工具更新，网络稳定也可以不设置

## Api Key准备
- 由于我在Fydeos 下使用的还是Google账户，所以申请google api key，不申请fydeos的，步骤见[官方说明文档](https://www.chromium.org/developers/how-tos/api-keys/)
- 大致步骤:
    + 加入chromium-dev开发群组
    + 登录`https://cloud.google.com/console`
    + 创建新项目，例如openfyde
    + 选择API服务 > 库，查找下面的所有API项目，逐一启用
        * Cloud Search API
        * Google Drive API
        * Safe Browsing API
        * Time Zone API
    + 进入API服务 > 凭据，根据提示先配置Oauth信息，然后创建凭据 > 创建 OAuth 客户端 ID，应用类型我选择了桌面应用（没有文档描述的Others）
    + 创建凭据 > 创建 API密钥，获得密钥
    + 保存客户端ID、客户端密钥、API密钥三个信息，按照下面方式存储到`~/.googleapikeys`文件中
        ```txt
        'google_api_key': 'your api key',
        'google_default_client_id': 'your client id',
        'google_default_client_secret': 'your client secret',
        ```

## 源代码获取

- 创建目录: `mkdir -p r120/openfyde`
- 配置Git，否则克隆代码会报错:
    ```bash
    git config --global user.name "sylens"
    git config --global user.email "silenseek14@gmail.com"
    ```
- 克隆源代码: R126版本的`openFyde/manifest.git`未更新，只能编译R120，`repo sync -j8`这里需要下载大量源代码内容，其中chromium最大，接近40G。如果网络问题导致下载失败，再次运行会重新下载，请确定梯子流量充足。
    ```bash
    repo init -u https://chromium.googlesource.com/chromiumos/manifest.git --repo-url https://chromium.googlesource.com/external/repo.git -b release-R120-15662.B

    git clone https://github.com/openFyde/manifest.git openfyde/manifest -b r120-dev
    ln -snfr openfyde/manifest .repo/local_manifests
    repo sync -j8
    cd openfyde/chromium
    gclient sync
    ```

## 进入chroot环境

- 回到`r120`这一层目录，运行`cros_sdk`，开始下载sdk并准备环境, 该过程时间较久

## 开始编译

- 由于当前版本archero无法使用，为避免报错，需要修改文件
- 进入chroot后，默认处在`/mnt/host/source/src/scripts`，需要修改的文件有:
    + `/mnt/host/source/src/overlays/overlay-fydetab_duo-openfyde/metadata/layout.conf`: 将第一行的`archero`、`tablet`、`ai-dev`删除
    + `/mnt/host/source/openfyde/overlays/overlay-fydetab_duo-openfyde/profiles/base/parent`: 将其中第二行的`archero:base`删除
- 注意: 编译的相关命令，都在chroot内执行
- 运行`setup_board --board=fydetab_duo-openfyde`
    + 看上去设置工作没有全部完成, 但是目标目录已生成
- 安装`capnproto`: `sudo emerge capnproto`
- 修改`/mnt/host/source/chromite/lib/dlc_allowlist.py`: 第14行修改为`DLC_FACTORY_INSTALL = (r"termina-dlc", r"sample-dlc",)`, 即加入`r"termina-dlc"`
- 继续运行`cros build-packages --jobs=4 '--board=fydetab_duo-openfyde' --no-withautotest --autosetgov --no-use-any-chrome`
- 运行完成后构建镜像`cros build-image --board=fydetab_duo-openfyde --noenable_rootfs_verification`

## 编译过程中的其他问题

- 编译异常中断: 遇到过一次异常中断，中断后直接进入环境再编译会有莫名其妙的问题，要从头开始需要删除`chroot`以及`out`文件夹，然后再次`cro_sdk`后从头再来。

- 虚拟硬盘设置: 经过[论坛大佬提示](https://community.fydeos.com/t/topic/48326/14)，还需要进行一项设置:
    + 打开`/etc/sysctl.conf`添加`vm.max_map_count=524288`
    + 运行`sudo sysctl -p`
