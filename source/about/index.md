---
layout: docs
seo_title: 关于我
bottom_meta: false
sidebar:
---

## 个人简介
- 自称: Silen Wang / Sylens Wong / 汪兴伦, 是一颗过得不是很开心的卤蛋!
- 学习/工作路线: 公卫 -> 分子流病 -> 生信 -> 新生抗原 -> 杂七杂八  emmmm... 越来越杂乱了...
- 兴趣: 打游戏, 打游戏, 还是TMD 打、游、戏! 动漫也看, 但是最近感兴趣番的越来越少了... 比起看番剧现在更喜欢看番剧吐槽-_-...


## 教育经历
2009.09-2014.06 广东药学院(现广东药科大学) 公共卫生学院 预防医学 本科
2014.09-2017.06 南方医科大学 公共卫生学院 流行病与卫生统计学 硕士


## 工作履历

- 2017.07-2018.08 天津某生信服务公司 生物信息分析工程师
    1. 人基因组二代测序数据（WES / WGS）质控及分析
        - 基于频率、基因区域、有害性软件预测、家系遗传模式进行变异结果筛选
        - 依据全基因组关联分析（GWAS）结果定位可能有害的变异 
        - 根据客户要求进行筛选方案定制并展示结果
    2. 特定疾病研究现状调研 
        - 对特定疾病的基本知识及研究现状进行调研学习，为开发相关分析模块作前置准备
    3. 数据交付流程维护与升级


- 2018.08-2019.07 杭州某基因科技公司 生信技术部 生物信息工程师
  - 数据处理
    1. cfDNA/FFPE/白细胞样品测序数据处理
    2. 二代测序分析流程构建 / 测试
    3. 二代测序分析用工具开发(Python / R / Julia / Shell)细胞样品测序数据拆分及质量控制
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

- 2019.07-2020.03 杭州某生物技术公司 生物信息工程师
  1. 基于NGS的cfDNA分析方案调研, 分析流程搭建
  2. 新生抗原相关数据统计分析
  3. 生物信息分析程序维护
  4. 用药指导报告数据库维护
- 2020.03-至今 杭州某生物技术公司 生信主管
  1. 基于NGS的HLA分型/定量模块开发
  2. 基于flask的生信分析工具后端开发
  3. 新生抗原筛选模块开发
  4. 生信分析主流程维护与升级


## 工作相关技能

- 脚本/编程语言/工具
    + 熟练: Python, R, Bash
    + 入门: 正则表达式, Graphviz, Mermaid, Go
    + 上手中: kivy
    + 已弃坑: Julia, LaTex

- 主力开发语言中有使用经验的模块
  - Python
    - 并行化: Multiprocess, Ray(上手中)
    - 二代测序常用文件处理: Pysam, PyVCF, CyVCF, Biopython
    - 爬虫: 
      + scrapy: 有数个小型爬虫项目经验, 能使用scrapy完成静态页面, 未加密动态页面的爬取和数据
      + 能使用splash, selenium配合scrapy完成动态页面数据的爬取
    - 数据处理/统计计算: Pandas, numpy, math
    - 绘图: Plotly
    - 数据交互展示: Dash
    - ~机器学习: skitlearn(上手中)~
    - 图像处理: pillow
    - web相关: flask / FastAPI
    - 办公自动化: docxtpl, openpyxl
  - R
    - base
    - 并行化: Parallel
    - 绘图: ggplot2, ggpubr, ggthemes
    - 单细胞数据分析: Seurat3

- 流程构建/部署:
    - Linux系统安装/配置管理: Redhat系发行版操作/维护, Debian系发行版安装/配置/维护, Arch系发行版安装/配置/维护
    - HPC集群管理系统: 有SUN Grid Eengine使用, 维护经验
    - Snakemake: 能熟练编写Snakefile, 能够使用Snakemake进行流程构建和管理, 有一定Snakemake的实际运用经验, 可以根据个人工作经验实现一些Snakemake暂时未实现的功能
    - ~~WDL~~: 由于WDL没有Python解释器, 暂时放弃了
    - Docker/Singularity: 有实际的容器部署经验, 能够通过docker-compose文件快速启动和部署应用. 能够通过撰写dockfile或直接进入容器内进行流程部署(展示项目中有示例), 可使用Singularity将Docker镜像转换为Singularity镜像并使用.
    - Conda/Mamba: 有较充足的实际使用经验, 能够使用conda快速部署分析用软件及流程.
    - Git: 有使用Git进行项目代码/文档管理的经验, 懂得基本的创建/推送/合并的方法, 使用过钩子特性运行一些自动化任务操作

- 掌握生物信息软件/模块(使用过及使用比较多的都有)
    - 二代测序数据拆分: bcl2fastq
    - 数据质控: fastp/fastQC/MultiQC
    - 序列比对: bwa/bowtie/STAR
    - 表达量计算/分析: DESeq/DEXseq/Salmon/kallisto/htseq/Hisat
    - 序列处理: seqtk/seqkit
    - SNP/INDEL检测: Samtools/GATK(Mutect2)/Vardict/Vacscan/Strelka
    - DNA结构变异/融合基因检测: CREAST/lumpy-sv/Smoove/SViCT/SvABA/Manta/Delly/GeneFuse/Factera/sv-tools
    - CNV检测: CNVkit/PyLOH
    - BAM/VCF文件操作: Picard/Samtools/Bcftools/Sambamba/Pysam/Biogo
    - BED文件操作: bedtools
    - 变异注释: ANNOVAR/snpEff/VEP
    - UMI处理相关: UMI-tools/fgbio
    - 克隆性: PyClone/SciClone/FastClone
    - 测序数据遗传一致性排查: plink/NGSCheckMate
    - 质谱鉴定: pFind/maxquant

- 其他工作经历
  - 专利稿撰写经历
  - 软著撰写经历
    
- 静态博客/文档撰写
  - hexo
  - mkdocs
  - gitbook


## 可展示项目

- [FUEX](https://github.com/SilenWang/FUEX): 融合基因解析软件, 从结构变异软件给的结果中解析出融合形式正确的结果
- [Risk Region](https://github.com/SilenWang/Risk_Region): 帮同学改写的(后续todo已鸽)
- [CSSN_Spider](https://github.com/SilenWang/CSSN_Spider): 帮同学做的简单爬虫, 爬了网站上所有的标准文档信息
- [DGIdb-Docker](https://github.com/SilenWang/dgidb-docker): DGIdb的docker构建文件