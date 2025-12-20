FROM mcr.microsoft.com/devcontainers/base:ubuntu
ARG PIXI_VERSION=v0.62.0

RUN curl -L -o /usr/local/bin/pixi -fsSL --compressed "https://github.com/prefix-dev/pixi/releases/download/${PIXI_VERSION}/pixi-$(uname -m)-unknown-linux-musl" \
  && chmod +x /usr/local/bin/pixi \
  && pixi info

USER vscode
WORKDIR /home/vscode
  
# 配置 git 用户信息
RUN git config --global user.name "Sylens" && \
  git config --global user.email "qiumin14@163.com" && \
  echo 'alias aider="aider --no-check-update --no-show-model-warnings --yes --no-auto-commits --model deepseek/deepseek-reasoner"' >> /home/vscode/.bashrc && \
  echo 'eval "$(pixi completion -s bash)"' >> /home/vscode/.bashrc