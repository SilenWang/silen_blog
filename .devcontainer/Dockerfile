# 使用 ms 提供的镜像，使用其他镜像会有用户问题（无vscode）
FROM mcr.microsoft.com/vscode/devcontainers/base:bookworm

# 安装常用工具
RUN curl -fsSL https://pixi.sh/install.sh | sh

# 配置 git 用户信息
RUN git config --global user.name "Sylens" \
  && git config --global user.email "qiumin14@163.com"