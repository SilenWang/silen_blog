---
title: mSINGS的安装和使用
categories: Bioinformatics
date: 2020-09-06 16:46:00
tags: ['微卫星不稳定性', 'MSI', '生物信息学', 'mSINGS', 'Python', 'samtools']
---

mSINGS是一个用来检测MSI的软件, 其优势似乎是可以用于tumor only的样品.
<!-- 摘要部分 -->
<!-- more -->

[mSINGS](https://bitbucket.org/uwlabmed/msings/src/master/)项目挂在bitbucket而不是github, 项目是一直都有在更新... 不过这个安装指引写的稍微硬核了一点... 然后它也没有一键安装之类的东西, 所以看着文档安装会有点绕...

这个项目的核心是用Python写的, 其核心是读取samtools生成的mpileup文件进行分析, 因此本质上的依赖应该是`Python3.6`和`samtools`就好了(里面的git应该是克隆项目用的吧?). 然后因为作者指定Python版本, 所以他才会建议使用Python虚拟环境(防止日后爆炸). 不过因为我不会用`virtualenv`, 所以用上手更简单的`miniconda`代替, 总的来说就是: 创建环境/安装必要依赖 -> 安装模块 -> 测试.

1. 创建环境/安装必要依赖

```bash
conda create -p /path/to/soft/mSINGS/conda python=3.6 git samtools
```

2. 安装模块

```bash
conda activate /path/to/soft/mSINGS/conda           # 进入环境准备安装模块本体
git clone https://bitbucket.org/uwlabmed/msings.git # 项目克隆, 这个按作者建议来
cd msings                                           # 进入目录
python setup.py install                             # 安装软件本体
```

3. 测试

作者在readme中提到了如何创建baseline(应该是根据已有数据确定哪些MSI位点需要扫), 但我只是测试, 所以直接用项目下准备好的文件进行(`doc/`目录下). 另外需要注意的是, 作者明确说了输入的bam需要比对到region不带`chr`字符串的参考基因组(GATK有提供的亚子), 所以我找了个fq数据从比对开始, 如果bam已经符合要求了可以跳过.

另外需要注意的是, 作者在运行脚本`run_msings.sh`, 指定了虚拟环境激活, 因为我没有按他的来, 所以脚本内的`# source msings-env/bin/activate`这一句需要注释掉.

```bash
set -e

# 不要做sort, 因为mSINGS的分析过程带了sort...
/path/to/apps/bwa-0.7.12/bwa mem \
    -R '@RG\tID:group1\tSM:TUMOR\tPL:illumina\tLB:lib1\tPU:unit1' \
    -M -t 16 \
    /path/to/Human_Genome_GRCh37_FASTA/human_g1k_v37.fasta \
    /path/to/soft/mSINGS/run_test/test.fq.R1.gz \
    /path/to/soft/mSINGS/run_test/test.fq.R2.gz \
| /usr/local/bin/samtools view \
    -Sb - \
> test.bam

echo "test.bam" > bam_list # 因为给的脚本指定一定要给写了bam路径的文件, 所以只能这么指定
sh /path/to/mSINGS/msings/scripts/run_msings.sh \
    bam_list \
    /path/to/mSINGS/msings/doc/mSINGS_TCGA.bed \
    /path/to/mSINGS/msings/doc/mSINGS_TCGA.baseline \
    /path/to/Human_Genome_GRCh37_FASTA/human_g1k_v37.fasta
```

最终的结果会在`Combined_MSI.txt`这个文件中, 这里面应该是将每个bam的结果放到了一起, 当然我只测试给一个bam, 所以合并的结果会不会有啥特殊结构未知...

我测试的这个样品目标MSI位点1166个, 27个检测到不稳定, 不稳定的2%, 总体是MSS.

```txt
Position        test
unstable_loci   27
passing_loci    1166
msing_score     0.0232
msi status      NEG
```

最后来看一下作者的运行脚本`run_msings.sh`:

```bash
#!/bin/bash
set -e
# source msings-env/bin/activate

#BAM_LIST is a file of absolute paths to each bam file
BAM_LIST=$1;
BEDFILE=$2;
MSI_BASELINE=$3;
REF_GENOME=$4;

#Check for required variables:
if [ -z "$BAM_LIST" ]; then echo "BAM_LIST is unset" && exit ; else echo "BAM_LIST is set to '$BAM_LIST'"; fi
if [ -z "$BEDFILE" ]; then echo "BEDFILE is unset" && exit ; else echo "BEDFILE is set to '$BEDFILE'"; fi
if [ -z "$MSI_BASELINE" ]; then echo "MSI_BASELINE is unset" && exit ; else echo "MSI_BASELINE is set to '$MSI_BASELINE'"; fi
if [ -z "$REF_GENOME" ]; then echo "REF_GENOME is unset" && exit ; else echo "REF_GENOME is set to '$REF_GENOME'"; fi


#"multiplier" is the number of standard deviations from the baseline that is required to call instability
multiplier=2.0 
#"msi_min_threshold" is the maximum fraction of unstable sites allowed to call a specimen MSI negative     
msi_min_threshold=0.2
#"msi_max_threshold" is the minimum fraction of unstable sites allowed to call a specimen MSI positive
msi_max_threshold=0.2

for BAM in `sed '/^$/d' $BAM_LIST`; do
    SAVEPATH=$(dirname $BAM)
    BAMNAME=$(basename $BAM)
    PFX=${BAMNAME%.*}

    mkdir -p $SAVEPATH/$PFX

    echo “Starting Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

    echo "sorting bam of $PFX" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    samtools sort -o $SAVEPATH/$PFX/$PFX.sorted.bam $BAM  && samtools index $SAVEPATH/$PFX/$PFX.sorted.bam

    echo "Making mpileups" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    samtools mpileup -f $REF_GENOME -d 100000 -A -E  -l $BEDFILE $SAVEPATH/$PFX/$PFX.sorted.bam | awk '{if($4 >= 6) print $0}' > $SAVEPATH/$PFX/$PFX.mpileup 
    
    echo "MSI Analyzer start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
    
    msi analyzer $SAVEPATH/$PFX/$PFX.mpileup $BEDFILE -o $SAVEPATH/$PFX/$PFX.msi.txt
    
    echo "MSI calls start" >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;
     
    msi count_msi_samples $MSI_BASELINE $SAVEPATH/$PFX -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/$PFX/$PFX.MSI_Analysis.txt

    echo “Completed Analysis of $PFX” >> $SAVEPATH/$PFX/msi_run_log.txt;
    date +"%D %H:%M" >> $SAVEPATH/$PFX/msi_run_log.txt;

done

echo "Creating summary analysis file for all samples" >> $SAVEPATH/msi_run_log.txt;
msi count_msi_samples $MSI_BASELINE $SAVEPATH -m $multiplier -t $msi_min_threshold $msi_max_threshold -o $SAVEPATH/Combined_MSI.txt
```

看着有点长, 其实内容还是比较简单的, 作者也给了一些基本的注释. 本质上就是将输入的bam一个一个去sort然后生成mpileup文件, 然后使用`msi`这个程序(之前`python setup.py install`装的东西)来进行分析和数据汇总.
