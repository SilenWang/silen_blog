# 使用 Pixi 官方基础镜像
FROM ghcr.io/prefix-dev/pixi:bookworm

# 安装常用工具
USER root
RUN apt-get update && \
    apt-get install -y sudo git curl bash vim ca-certificates less && \
    rm -rf /var/lib/apt/lists/*

# 创建 vscode 用户并允许无密码 sudo
RUN useradd -ms /bin/bash vscode && \
    echo "vscode ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# 设置工作目录并切换用户
WORKDIR /workspace
USER vscode
