# 使用 ms 官方基础镜像，含有vscode用户，避免devpod启动时有问题
FROM mcr.microsoft.com/devcontainers/base:bookworm

# 为vscode安装pixi
USER vscode

# 安装常用工具
RUN curl -fsSL https://pixi.sh/install.sh | sh

# 配置 git 用户信息
RUN git config --global user.name "Sylens" && \
  git config --global user.email "qiumin14@163.com" && \
  echo 'alias aider="aider --no-check-update --no-show-model-warnings --yes --no-auto-commits --model deepseek/deepseek-reasoner"' >> ~/.bashrc

ENV PATH="/home/vscode/.pixi/bin:${PATH}"