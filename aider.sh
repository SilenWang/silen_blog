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
        echo "开始文稿校对"
        /home/sylens/.pixi/envs/pip/bin/aider --no-show-model-warnings --yes --no-auto-commits \
            --no-show-model-warnings --yes --no-auto-commits \
            --message "我这里有一篇文稿的中文版本，请根据中文版本，生成英文版本的文稿，英文版请存放到 ${DIR_B}/${filename} 中" \
            $file
    fi
done