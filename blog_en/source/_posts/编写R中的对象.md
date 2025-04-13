---
title: Writing Objects in R
categories: Bioinfomatic
date: 2025-04-13 15:57:55
tags: ['rlang', 'object']
---

Object-oriented programming is a common paradigm in programming. In my actual work, using object-oriented programming mainly serves two purposes: reducing duplicate code through inheritance features, and encapsulating frequently used data into objects to avoid excessive, repetitive, and nested parameter passing.

<!-- Excerpt -->
<!-- more -->

## Different Object Systems in R

Yes... R has different object systems. Besides the previously encountered S3 (like ggplot) and S4 (adopted by many single-cell packages), there's R6 which we'll focus on this time.

S3 is the most basic and widely used object-oriented system in R, implementing method dispatch through **generic functions**. S3 objects essentially feel like extensions of lists, with relatively simple inheritance implemented via generic functions. This differs significantly from the class systems I'm more familiar with in Python.

S4 provides a more formal object-oriented programming framework suitable for complex applications requiring strict class checking. It's more complex than S3 with more comprehensive inheritance mechanisms, but its methods are still implemented through generic functions. While still different from Python's approach, with inheritance and methods available, and the ability to conveniently use the `%>%` pipe operator, S4 is somewhat usable.

Originally I wanted to base my development on S4, but unexpectedly, S4 and `box` have [compatibility issues](https://github.com/klmr/box/issues/284). The author of `box` stated in 2022 that it wouldn't be compatible with S4...

So finally, I used the `R6` class, which is a third-party package (yet another third-party one...). The objects it provides are very similar to Python objects - with inheritance, methods, and support for `box::use`. The only drawback is I haven't found a way to integrate it with `%>%`.

## Using R6 Objects

Here's a concrete example of using R6, which indeed closely resembles Python's approach, though R6 doesn't require defining an initialization function.

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

# Using the object
person <- Person$new("John", 30) # Instantiation
person$greet() # Calling the method
```
