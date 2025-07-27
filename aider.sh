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
            --message "我这里有一篇文稿的中英双语版本，请校对两者的内容是否一致，如果不一致，按照中文的描述来修正英文文稿中的描述，然后将修改保存到原文件中" \
            $file ${DIR_B}/${filename}
    fi
done