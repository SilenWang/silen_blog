---
title: python中的itertools
categories: Coding
date: 2019-02-13 00:37:27
tags: ["Python", "itertools"]
---

I've never used many useful functions when writing scripts before, so I decided to learn a bit about `itertools` and make some notes.

This time, I was dealing with large files, but most of the time, we handle them line by line. So, I thought about using `itertools` to create generators. This way, I might be able to split the generator into smaller parts and apply it to parallel operations.

The three functions from `itertools` that I used this time are:

1. groupby
2. tee
3. chain

- `groupby` is used to cluster items in an iterable object. For example, when dealing with UMI grouping, you need to compare reads that are close in position and group them together for further UMI analysis.

```python
# Here's a function that takes a single object from iter_obj as input and returns the value to be used for comparison
grp_generator = groupby(iter_obj, lambda x: x[0])
# Note that groupby only groups adjacent objects, so you might need to sort them yourself.
```

- `tee` is simple. If the iterable object is a generator, since generators can only be used once, you can use the `tee` function to create a copy of it for multiple uses.

```python
iterable_fork_1, iterable_fork_2 = tee(iterable, 2) # The number 2 is optional and defaults to 2.
```

- `chain` concatenates multiple iterable objects into one. I used it to flatten nested results:

```python
flat_iter = chain.from_iterable(ori_iter)
# [[1,2],3] -> [1,2,3] Of course, the actual result is a generator, so you need to use list(flat_iter) to get the complete result.
```
