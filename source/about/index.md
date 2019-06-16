---
title: 关于我
date: 2019-03-11 00:03:29
---

## 个人简介
公卫 -> 流病 -> 生信 -> 算法 emmmm... 总觉的我的路线是不是有点歪了...


## 教育经历
2009.09-2014.06 广东药学院 公共卫生学院 预防医学 本科
2014.09-2017.06 南方医科大学 公共卫生学院 流行病与卫生统计学 硕士


## 工作履历


- 2017.07-2018.08 天津诺禾致源科技有限公司 人类疾病业务线 生物信息分析工程师
    1. 人基因组二代测序数据（WES / WGS）质控及分析
        - 基于频率、基因区域、有害性软件预测、家系遗传模式进行变异结果筛选
        - 依据全基因组关联分析（GWAS）结果定位可能有害的变异 
        - 根据客户要求进行筛选方案定制并展示结果
    2. 特定疾病研究现状调研 
        - 对特定疾病的基本知识及研究现状进行调研学习，为开发相关分析模块作前置准备
    3. 数据交付流程维护与升级


- **2018.08-至今** 杭州瑞普基因科技公司 生信技术部 生物信息工程师
  - 数据处理
    1. cfDNA/FFPE/白细胞样品测序数据处理
    2. 二代测序分析流程构建 / 测试
    3. 二代测序分析用更具开发(Python / R / Julia / Shell)细胞样品测序数据拆分及质量控制
    4. 原始测序数据(fastq)比对, 及质量控制
    5. UMI数据 / Molecular Barcode数据处理
  - 分析流程相关
    1. 对现有生物信息分析流程Debug
    2. 编写新的分析模块并整合进流程框架
    3. 分析模块的测试及Debug
    4. 分析用软件部署及依赖处理
    5. 分析用软件性能测试及分析结果比较
    6. 流程文档的编写与维护
    7. 单细胞测序分析流程的搭建及测试
  - 算法及工具开发:
    1. UMI去重算法比较及复现
    2. 质控图形绘制工具开发
    3. 基于机器学习的变异结果筛选工具的开发
    4. 已有分析工具的效能优化(并行化)


## 掌握技能


- 脚本/编程语言/工具
    + 熟练: Python, R, Bash
    + 入门: 正则表达式, Graphviz, Go
    + 上手中: Julia, LaTex


- 主力开发语言中有使用经验的模块
  - Python
    - 并行化: Multiprocess, Ray(上手中)
    - 二代测序常用文件处理: Pysam
    - 常用格式处理/常用数据库爬取: Biopython
    - 数据处理/统计计算: Pandas
    - 绘图: Plotly
    - 机器学习: skitlearn(上手中)
    - 图像处理: pillow
  - R
    - base
    - 并行化: Parallel
    - 绘图: ggplot2, ggpubr, ggthemes
    - 单细胞数据分析: Seurat3


- 流程构建/部署:
    - Linux系统安装/配置管理(Redhet系列, Debian系列, Arch系列)
    - Snakemake: 能熟练编写Snakefile, 能够使用Snakemake进行流程构建和管理, 有一定Snakemake的实际运用经验, 可以根据个人工作经验实现一些Snakemake暂时未实现的功能
    - ~~WDL~~: 由于WDL没有Python解释器, 暂时放弃了
    - Docker: 有实际的容器部署经验, 能够通过docker-compose文件快速启动和部署应用. 能够通过撰写dockfile或直接进入容器内进行流程部署.
    - Conda: 有较充足的实际使用经验, 能够使用conda快速部署分析用软件及流程.
    - Git: 有使用Git进行项目代码/文档管理的经验, 懂得基本的创建/推送/合并的方法


- 掌握生物信息软件/模块(使用过及使用比较多的都有)
    - 二代测序数据拆分: bcl2fastq
    - 数据质控: fastp/fastQC
    - 序列比对: bwa/bowtie
    - SNP/INDEL检测: Samtools/GATK(Mutect2)/Vardict
    - DNA结构变异/融合基因检测: CREAST/lumpy-sv/Smoove/SViCT/SvABA/Manta/Delly/GeneFuse/Factera/sv-tools
    - BAM/VCF文件操作: Picard/Samtools/Bcftools/Sambamba/Pysam/Biogo
    - BED文件操作: bedtools
    - 变异注释: ANNOVAR/snpEff
    - UMI处理相关: UMI-tools/fgbio
    
    
- 静态博客/文档撰写
  - hexo
  - mkdocs
  - gitbook


## 个人项目

- [FUEX](https://github.com/SilenWang/FUEX): Fusion extracter for SV detection tools