---
title: R Small Tips Collection
categories: Coding
date: 2018-09-07 12:01:44
tags: ['R']
---

Working with R often involves encountering various issues. Here are some collected tips.
<!-- more -->

# Data Frame Operations

## Extracting a Single Column/Row Without Returning a Vector

- When using `df[]` to subset a data frame, if you extract a single column or row, or delete all but one column or row, the resulting object will be a vector instead of a data frame. You need to specify the parameter to ensure it returns a data frame.
- This method also applies to matrices.

```r
df[c(1, 2, 3), -2, drop = F]
```

## Merging Data Using row.names as Key

- Many Chinese materials describe merging data by using one column as the key for row-wise merging. However, you can also use `row.names` as the key, but after merging, `row.names` will become a column of data. You need to specify `row.names` again.

```r
# Two-step process
merge_data <- merge(x, y, all.x = T, by = 0)
merge_data <- data.frame(merge_data), row.names = 1)
# One-step merge
merge_data <- data.frame(merge(x, y, all.x = T, by = 0), row.names = 1)
```

# Other

## R's Python Dict Alternative

- Many online resources use lists combined with the `match` function to implement dictionary-like functionality. However, in R, named vectors can be used to achieve key-value mapping, which is very convenient. However, it's unclear whether this approach is efficient.

```r
a <- c(1, 2, 3)
names(a) <- c("a", "b", "c")
tmp <- a["a"]
print(tmp)
```
```