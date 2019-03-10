---
title: snakemake的使用
categories: Bioinfomatic
date: 2018-10-30 12:40:03
tags: ['python', 'snakemake']
---

写流程是生物信息分析中经常需要涉及的工作, 一个成熟完善的好工具可以使工作的效率大大提高.
<!-- more -->

最简单的流程编写方式可以是将所有需要运行的命令都写成一个shell脚本, 并使用`&&`或其他方式指定任务的先后顺序, 然后挂到后台就好.

这么做虽然容易, 但是其重复性和抗错误性实在堪忧. 因此也就有了各种各样的流程管理工具: 我过去单位用的`sjm`, 10x Genomics自己开发的`Martian`, Broad Institute的`WDL`, 我最近在用的`snakemake`等等.

## snakemake简介
`snakemake`是德国杜伊斯堡-埃森大学的Köster lab开发的生物信息分析流程管理工具, 其全部代码完全由python实现, 因此可以很好的和python程序整合在一起(当然其也有一套完整的命令行程序).

`snakemake`的特点包括但不限于:

- 简单易用
- 可直接使用python代码
- 内置分布式运算集群支持(SGE, LSF, etc.)
- 内置Conda环境部署支持
- 内置容器支持(注意这个容器不是docker)
- 自动处理依赖并进行断点重启
- 内置测试功能(benchmark)
- 与`make`完全一样的执行逻辑, 有相关基础的人可快速上手

## snakemake简单示例

snakemake文件中以`rule`为最基本单位, 一个`rule`就算是流程中的一个步骤. 在`rule`下可定义输入文件, 输出文件以及要执行的命令. 由于snakemake基于python, 所以其在语法格式上也与python一致, 也是以缩进来定义控制范围, 比如写一个简单的`rule`:

```python
rule samtools_sort:
    input: 
        "test.bam"
    output:
        "test.srt.bam"
    shell:
        "samtools sort {input} > {output}"
```

这个名称为"samtools_sort"的`rule`执行的步骤是将输入的`test.bam`通过`samtools`的`sort`命令进行排序, 然后将输出结果指定为"test.srt.bam".

最后的`shell:`代表这里执行的是shell命令, 命令以字符串的形式给出, snakemake还支持直接执行python代码, 此时的关键词需要更改为`run:`. 此外也支持调用R代码以及其他外部脚本, 详情参考snakemake的[官方手册]().

## 使用snakemake进行样品批量处理

了解了最基本的`rule`写法之后就可以写出一串`rule`, 然后把他们串成流程了, 不过虽然snakemake简单易学, 要写出一个能处理一堆样品的完整流程还是需要了解些东西的.

### snakemake中的特殊字符串格式化方式

从之前的简例中可以看出, snakefile中执行的命令, 输入及输出的文件都以字符串的形式给出. 而在这些字符串中可以使用类似`str.format()`的格式化字符串方式减少硬编码, 增加应用性.

比如前面的`"samtools sort {input} > {output}"`中, `{input}`会被替换成**该rule内定义的**`input`后面字符串(如果input不只一个, 则会以空格分隔全部列出). 类似的, 我们可以调用任意`rule`内定义过的内容来填充文本模板, 以生成最终需要运行的命令或者字符串.

但是有一点需要注意, 对于习惯了用`format()`的人来说可能会写类似下面的表达:

```python
path="path/to/input"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {path}/{input} > {output}".format(path=path)
```

即将python中的字符串格式化和snakamake的字符串格式化混用. ~~这种写法snakemake是没法正确处理的, 一个字符串要么只能使用内置的格式化方式, 要么只能用python原有的格式化方式, 上面那种写法, 会按照python原有的方式去处理, 即根据`format()`内给的信息去填充字符串, 没有指定的`input`和`output`会作留空处理.~~ snakemake是可以处理这种情况的, 但是写法上要进行变更:

```python
path="path/to/input"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {path}/{{input}} > {{output}}".format(path=path)
```

用python比较多的话, 一定了解`"\{\{\}\}"`代表不对字符串内的`{}`进行转义, 而当作一般符号. 结合snakemake的执行逻辑, 可知snakemake会自行处理一次纯字符串的输入, 而这一切发生在python部分执行完之后.

### snakemake的流程执行逻辑

snakemake并没有准备专门的顺序或依赖的语法关键词, 流程的串联完全依靠输入和输出文件的依赖关系自动完成. 由于snakemake默认会去解析并完成第一个`rule`, 因此官方文档推荐的使用方式是创建一个名称为"all"的`rule`, 然后在这个`rule`的`input`中(注意不是`output`)指定所有的最终文件. 然后程序会发现结果文件不存在, 便在文件中寻找到某条`rule`的`output`符合文件模式, 然后继续解析这个`rule`需要要什么作为`input`, 如果`input`(即执行依赖)已满足, 则会开始运行该规则下的命令/脚本, 如果不满足, 则会继续向前查找, 直到找到源头或者找不到源头而报错为止, 例如:

```python
rule all:
    "test.flt.srt.bam"

rule samtools_flt:
    input: 
        "test.bam"
    output:
        "test.flt.bam"
    shell:
        "samtools sort {input} > {output}"

rule samtools_sort:
    input: 
        rules.samtools_flt.output
    output:
        "test.flt.srt.bam"
    shell:
        "samtools sort {input} > {output}"
```

在上面的snakefile中, 流程运行的逻辑是:

1. 最终目标是`test.flt.srt.bam`, 搜寻可以产生这个文件的`rule`
2. 发现`samtools_sort`可产生需要的文件, 而这个规则需要输入`samtools_flt`规则的结果, 因此`samtools_flt`先执行
3. `samtools_flt`需要输入`test.bam`, 查询该文件是否存在, 如果存在则按照`samtools_flt > samtools_sort > all`的顺序执行流程, 如果找不到该文件, 则会提示依赖不满足, 流程无法运行.

由于这种回溯式的流程执行逻辑, 在流程中某一个`rule`失败后, 非常容易实现流程断点重启: 即回溯式的检查每个步骤的依赖是否已满足, 如果满足, 则从当前步骤执行, 如果不满足, 则继续向上回溯.

另外, snakemake在和执行流程的过程中, 本身会检查输出文件的更新情况以及文件完整性以确定规则是否顺利运行. 如果某个规则出错, 改规则下的结果文件默认会被删除, 所以除非snakemake本身的进程被异常终止, 否则不太会出现某个步骤未完成但是流程继续向后跑的情况.

### snakemake内建的实用函数

由于snakefile内是可以直接写python代码的, 所以对于用惯了python的人, 可以比较轻松的使用生成式之类的东西快速生成批量处理时需要的信息. 而对于没有太多python基础的人, snakemake也准备了一些内建的函数来快速实现类似的功能. 比如:

- expand: 
    + 用来生成字符串列表的函数, 相当于组合了列表生成试和格式化字符串, 默认情况下会将多个替换变量两两组合, 但也可以设置为类似`zip()`的模式, 即依次组对, 可以满足不同情况下的需要
    + `samples = expand("{sample_id}.{fq}.fasta", sample_id=['sam1', 'sam2'], fq=['fq1', 'fq2'])`

### 多样品/多文件的批量处理

批量处理涉及到snakemake中的一个比较重要的内置对象----`wildcards`. 我虽不是太清楚这个东西应该怎么翻译...大致上来说, 这个对象用来进行种带捕获效果的模式匹配(正则). 比如官方文档中给出的例子:

```python
rule complex_conversion:
    input:
        "{dataset}/inputfile"
    output:
        "{dataset}/file.{group}.txt"
    shell:
        "somecommand --group {wildcards.group} < {input} > {output}"
```

这里规则的输入输出中的`{dataset}`以及`{group}`就是两个`wildcard`, 而每一个`wildcard`都可以匹配近乎任意字符(对应正则表达式为`.+`)如果某个规则需要一个`file_path/file.001.txt`的文件. snakemake在进行查找的时候就会发现这个文件的命名符合`complex_conversion`下的`output`规定的模式, 因此snakemake认为`file_path/file.001.txt`由该规则产生, 并将`dataset`设定为"file_path", `group`设定为"001", 然后把这两个变量向上传递到`input`, 继续进行依赖解析.

如此依赖, 当请求的结果文件是按照一定规则去命名的(通常来说会有样品/项目ID之类的信息), 就可以通过`wildcard`来匹配每个样品/项目需要进行哪些步骤, 然后自动生成完整的流程路线并开始执行.

## snakemake & conda

为了保证已完成流程的可重复性, snakemake内置了使用conda来进行流程依赖部署的功能.
不过这个功能和我最开始理解的不太一样.

## snakemake & 容器

除了conda, snakamake还支持使用容器来进行依赖部署. 只不过这个容器并非最火热的docker, 而是singularity. 我也是用了snakemake才知道原来docker只是容器计数中的一种实现, 只是它特别的广为人知. snakemake支持的singularity与docker相比的特点是:

- 不需要root权限
- 更低的资源占用

该项目从3.0版本开始从c++迁移到了go, 我虽然克服了网络压力装上了, 但应该是由于发行版的问题不能正常使用, 所以暂时没有实际的使用经验.

## 问题解决方案收集

本人在使用snakemake过程中时不时会碰到一些大小问题, 由于这个软件的中文使用者似乎并不多, 英文材料虽然完善但是由于我相关的储备不足, 并不能很快从中找到解决方案, 因此将我碰到的问题&解决方案收集在此, 部分问题的描述可能并不专业, 以后有机会慢慢更正