set -e 

# DIR_A=blog_cn/source/_posts
# # DIR_B=blog_en/source/_posts

# for file in `ls $DIR_A/*`; do
#     # 获取文件名（去掉路径部分）
#     filename=$(basename "$file")

#     # 检查文件是否存在于目录B
#     echo "开始对 ${filename} 中的内容进行关键词整理"
#     /home/sylens/.pixi/envs/pip/bin/aider \
#         --no-show-model-warnings --yes --no-auto-commits \
#         --message "这里是一篇 hexo 框架的博客文稿，请根据文稿的内容，提炼这个文稿的关键词，然后替换现在文稿中规定的关键词。请注意关键词的排列顺序，中文的关键词统一放前面，英文的关键词则放后面" \
#         $file ${DIR_A}/${filename}
# done


filename="blog_cn/source/_posts/我是什么天煞孤星体质么.md"

# 检查文件是否存在于目录B
echo "开始对 ${filename} 中的内容进行关键词整理"
/home/sylens/.pixi/envs/pip/bin/aider \
    --no-show-model-warnings --yes --no-auto-commits \
    --message "这里是一篇 hexo 框架的博客文稿，请根据文稿的内容，提炼这个文稿的关键词，然后替换现在文稿中规定的关键词。请注意关键词的排列顺序，中文的关键词统一放前面，英文的关键词则放后面" \
    ${filename}