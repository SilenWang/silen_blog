---
title: 编写R中的对象
categories: Bioinfomatic
date: 2025-04-13 15:57:55
tags: ['rlang', 'object']
---

面向对象是编程时一种常用的范式，在我的实际工作中，使用面向对象主要是为了通过继承特性减少重复代码，和将常用数据封装到对象内，避免过多、重复、嵌套的传递参数。

<!-- 摘要部分 -->
<!-- more -->

## R中不同的对象系统

是的... R中有不同的对象系统，除了之前接触的S3（如ggplot）和S4（被很多单细胞包采用），还有本次着重说明的R6。

S3是R语言中最基础、使用最广泛的面向对象系统，采用**泛型函数**实现方法分派。S3的对象，感觉本质上是list的延伸，它的继承比较简单，方法是通过泛型函数来实现的，跟我更熟悉的python中的类型差别较大。

S4提供更正式的面向对象编程框架，适合需要严格类型检查的复杂应用，它比S3要复杂，且有更全面的继承机制，但是其方法同样是通过泛型函数来实现的。虽然还是跟python中的不太一样，但是有继承，有方法，同时使用S4也能比较方便的使用`%>%`管道，勉强能用。

本来我想基于S4来完成我的开发，但是，没想到S4和`box`居然还[有兼容性问题](https://github.com/klmr/box/issues/284)，`box`的作者2022年就说明，不会兼容S4...

于是最后，我使用了`R6`类型，它是第三方的包（又是第三方...），其提供的对象，非常类似python的对象，有继承、有方法，支持`box::use`，美中不足，还没有找到接入`%>%`的方法。

## 使用R6对象

R6的具体使用实例如下，方式确实非常类似python，不过R6不需要定义初始化函数。

```r
box::use(
    R6[R6Class]
)

Person <- R6Class(
    "Person",
    public = list(
        name = NULL,
        age = NULL,
        initialize = function(name, age) {
            self$name <- name
            self$age <- age
        },
        greet = function() {
            cat("Hello, my name is", self$name, "and I am", self$age, "years old.\n")
        }
    )
)

# 使用对象
person <- Person$new("John", 30) # 实例化
person$greet() # 调用方法
```
