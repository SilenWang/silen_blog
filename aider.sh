set -e 

DIR_A=blog_cn/source/_posts
DIR_B=blog_en/source/_posts

for file in `ls $DIR_A/*`; do
    # 获取文件名（去掉路径部分）
    filename=$(basename "$file")

    # 检查文件是否存在于目录B
    if [[ -e "${DIR_B}/${filename}" ]]; then
        echo "文件 ${filename} 存在于目录${DIR_B}，跳过"
    else
        echo "文件 ${filename} 不存在于目录${DIR_B}，启动aider翻译..."
        aider \
            --no-show-model-warnings --yes --no-auto-commits \
            --model ollama_chat/hhao/qwen2.5-coder-tools:7b \
            --message "请翻译 ${file} 的内容为英文，然后将翻译后的文本保存到 ${DIR_B}/${filename}" \
            $file
    fi
done