---
title: pandas Common Operations Record
categories: Coding
date: 2020-07-26 00:36:49
tags: ['Pandas', 'Python']
---

<!-- Abstract part -->
<!-- more -->
As my current work requires a lot of data statistics and processing, and considering the convenience in error handling and third-party extensibility, Python is much more convenient than R. Therefore, I have started moving my data processing from R to Pandas. In this post, I will record some common operations used in actual work and note down any issues encountered.

## 1. apply Related

### 1.1. Basic Usage

Similar to Pandas and R, there is only one `apply` function, no `apply` family. The `apply` function essentially loops through the table row by row. Using `apply` instead of directly `itterrow` is mainly for parallel considerations. There are many third-party modules that can parallelize `apply` operations, so initially writing code in this way saves the trouble of rewriting functions later when parallelizing.

Of course, another point is to separate the main processing logic from detailed execution content, which to some extent improves code readability. However, personally, it's more of a habit issue... If the person reading the code never writes like this... Then using `apply` actually increases the difficulty of reading the code for others...

Currently, my approach is generally to write the actual calculation part as a general function and then use an anonymous function to pass parameters when calling `apply`.

```python
def plus(x, y):
    return x + y

data['applyed'] = data['2apply'].apply(lambda x: plus(x['x'], x['y']))
```

Of course, if you need to perform many row-wise operations at the same time, it's better to write a separate `func_apply()` function to prevent repeated traversal of the data.

### 1.2. Returning Multiple Columns

The conventional usage of `apply` is to return a numerical value to form one column of data. However, in actual applications, I often need to return multiple columns. Therefore, I found a solution online:

```python
def plus(x, y):
    return x + y, x - y


data['x+y'], data['x-y'] = zip(*data['2apply'].apply(lambda x: plus(x['x'], x['y'])))
```

## 2. Groupby Related

### 2.1. Merging Multiple Rows into a Single Row of Data

This operation is mainly used to aggregate multiple rows of data into a single row result. I remember that R has a specific way to handle this, and Pandas should have it too... However, I couldn't find it at the time, so I came up with my own solution. This method essentially retrieves data subsets based on the unit to be processed, performs statistics on these subsets, and generates new rows of data to fill in a new table. This approach currently seems inefficient, and because it uses a for loop, parallelization is not as simple as `apply`.

```python
result_list = []

for k, sub_data in data.groupby(by=['Grp1', 'Grp2']):
    grp1, grp2 = k
    result_list.append({
        'Stat1': sum(sub_data['Var1']),
        'Stat2': min(sub_data['Var1']),
    })
result = pd.DataFrame(result_list)
```

## 3. Missing Value Handling Related

### 3.1. Missing Values on Read

You can replace multiple values with missing ones using the `na_values` parameter:

```python
df = pd.read_table(FILE, sep="\t", header=0, na_values=['.', ''])
```

### 3.2. Missing Values on Output

When outputting, you can specify what to fill in for missing values through parameters:

```python
df.to_csv(FILE, sep="\t", header=0, na_values=['.', ''])
```

### 3.3. Determining if a Value is Missing

In Pandas, if a missing value is generated, it uses the `nan` object from NumPy, which is different from R's `NA`, `NULL`, and Python's built-in `None`. The type of `nan` is actually `float`... Therefore, using logical judgment like `A == nan` will not get `True` results... I haven't found a reason for this. However, to determine if a value is missing, it's best to use the `isnull` method provided by Pandas. Other modules like `math`, `numpy` have built-in missing value judgment functions, but these functions have a fatal flaw... They only accept parameters of type `float`... So if I use them, I still need to add an additional type judgment layer, which is quite troublesome. Therefore, using the `isnull` method provided by Pandas is the most convenient.

## 4. Most Common Output Parameters

Because I habitually output TSV files, when using Pandas to output files, I usually write it like this:

```python
df.to_csv(FILE_PATH, sep="\t", index=0, na_rep=".")
```

It's worth noting that if the result table is a pivot table, there is a high probability of needing to output the `index`. In such cases, `index=0` should be removed.

## 5. Other Possible Operations

### 5.1 Column Type Conversion

```python
df.astype({'col1': 'int32'})
```
```