---
title: R中的逻辑回归
categories: Others
date: 2019-03-16 14:07:19
tags: "R"
---

想不到4年之后, 又要弄回归了...这一次用的是R而不是SPSS了
<!-- 摘要部分 -->
<!-- more -->

## 简介
逻辑(logistic)回归是广义线性模型的一种, 其特点在于因变量是分类变量, 通过将所有分类中的一个设为对照, 将其他的分类选项发生的的可能都与对照做比较的方式, 将分类变量转变成连续变量, 然后进行模型构建. 在逻辑回归中, 变量的系数取`exp(x)`后即为这个变量对应的OR值, 实际意义为, 该变量每增加一个单位, 接过被分类为
1的可能性将提高`exp(x)`.

## 模型构建
本次我使用的逻辑回归因变量是二分类的, 所以使用基本包的`glm()`函数就够了, 具体模型构建使用:

- `family=binomial(link='logit')`指定调用二分类逻辑回归模型

```r
model <- glm(y ~ x1 + x2, family=binomial(link='logit'), data=fit_data)
summary(model) # 显示模型信息
```

在回归的时候, 很难保证第一次回归就能得到最好的结果. 在纳入的协变量中, 总会有系数过小的(也就是实际没太大作用), 或者统计检验不显著的(再抽样一次可能就没作用的), 因此我们需要变更变量多次进行回归.

在R中`step()`可以帮助我们自动进行回归, 并将对模型贡献小和不显著的变量自动去除:

```r
model.step <- step(model)
summary(model.step) 
```

当然, 这仅仅是个自动化的参考, 最后模型怎么确定好还是要根据实际情况由人来进行具体判断.

## 模型预测情况估计

模型构建好后, 即可以用`predict()`对测试数据集做预测, 以对模型进行评价, 然后使用一些第三方包通过ROC曲线确定拟合情况, 并对最终的分类阈值做选择

```r
test_data$prob <- predict(model.step, test_data, type = "response") # response的话接过是0~1的概率, 不指定默认是个常数

# 使用第三方的pROC
library(pROC)

obj_roc <- roc(test_data$Real, test_data$prob) # Real是分类的正确答案, prob是给的预测概率
plot(obj_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), # 做图
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
# pROC给的ROC曲线图y轴是灵敏度(真阳性), x轴则是特异度(真阴性), 而非更常见的假阳性率, 所以x轴是1->0而非0-1
```

## 小Tips

生成的模型对象本身也会包含详细的模型信息, 但是直接从模型上取信息与`summary()`看到的并不一致, 如系数的统计表, `summary()`给出的表格包含了系数数值, 系数表准误, 统计检验P值等详细信息, 但是直接从model取的话就只有系数数值.

想要取`summary()`信息的话, 要在其生成的对象上取, 或者先将其赋值后再取:

```r
# 先赋值
model_info <- summary(model)
model_info$coefficients
# 直接取
summary(model)$coefficients
```