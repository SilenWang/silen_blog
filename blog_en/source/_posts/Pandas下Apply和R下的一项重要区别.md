---
title: A Key Difference Between Pandas in Python and R
categories: Coding
date: 2023-01-15 00:30:41
tags: ['Pandas', 'R', 'Python']
---

I currently rarely use R for data processing or cleaning because my daily work involves a lot of string extraction/processing, which is difficult to do in R. Additionally, the error tracking in R is very demanding on code proficiency and experience, making it almost unbearable for me to write and maintain R code (I tried at a previous company...). However, recently one of my colleagues asked me to use R to complete such tasks because he only knows R. In the process of copying code, I discovered another reason why I don't need to use R for this type of work...

<!-- Abstract part -->
<!-- more -->

The content that needs to be processed is essentially taking a subset of the data frame and then generating new columns based on the row content of the subset, before concatenating these subsets. I have done this countless times with Pandas... Therefore, I gradually mastered the use of `apply`, `group`, `pivot_table`, and `melt` functions in Pandas, and there is nothing that these functions cannot solve for me, except for my poor coding ability, which makes the code inefficient...

However, performing the same operation in R is much more difficult... After all, it's two different languages. Trying to find an equivalent or similar function in R to handle data as I do with Pandas is very challenging, and even when there are equivalents like `apply`, they still cause problems that I cannot solve, leading me to revert back to using a `for` loop...

The processing to be done involves selecting a subset `sub` from the dataset `data`. This `sub` actually contains two parts: `A` and `B`. The task is to perform row-wise judgment, returning part `A` if the content of `A` is not missing, otherwise returning part `B`, with each row's return value forming a new data frame.

This can be easily achieved in Pandas by applying a function to the subset, then writing a function to judge the row and return the result. After `apply` completes, you directly get the desired data frame. However, there is a huge problem here: the data types of elements within a vector must be the same. When I collected numerical values into the same vector, non-string elements were automatically converted to strings, so although it could run, the newly generated four columns are all strings.

To try and solve this issue, I replaced vectors with `list` in the function, resulting in an intermediate version. Then I thought of another problem: `t()` transposes data as a matrix because R's columns are individual vectors. If you don't do it this way, transposing will result in one column having multiple data types... Therefore, I searched online and copied a piece of code to avoid using `t()`, resulting in the following version:

```R
# Function to return the required content
mkSel <- function(rec){
  mkA <- rec[['MarkerA']]
  mkB <- rec[['MarkerB']]
  if (is.na(mkB)) {
    row <- list(mkA, rec[['ColA1']], rec[['ColA2']], rec[['ColA3']])
  } else {
    row <- list(mkB, rec[['ColB1']], rec[['ColB2']], rec[['ColB3']])
  }
  names(row) <- c('Marker','Tag1', 'Tag2','Tag3')
  return(row)
}

sub[,c('Marker','Tag1', 'Tag2','Tag3')] <- do.call(rbind.data.frame, apply(sub, 1, mkSel, simplify=T))
```

However, after testing, I discovered the biggest problem... `apply` function... it simply cannot return a data frame; it can only return a matrix, and matrices are inherently of the same type... Even the entire family of functions does not have one that can meet my needs. So no matter how hard I tried, using `apply` to achieve my goal was impossible.

In the end, I reverted back to the humble `for` loop to solve the problem... Generate an empty data frame, specify the types for each column, read each row of the `sub` data frame to generate a new data frame, and finally merge the two data frames according to their indices...

So in summary, R's `apply` was not designed with this usage in mind; it can only use basic functions... Through a more cumbersome (needing to write several lines at a time, while also seeming less efficient and less reusable) way of achieving what I do with Pandas...

Of course, objectively speaking, my proficiency in R is far inferior to that of Pandas. I haven't used many useful R packages (like `dplyr`), so maybe R can complete these tasks well... But... I couldn't find any relevant materials?
