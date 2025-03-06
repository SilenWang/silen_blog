---
title: Pandas Operation Notes
categories: Script
date: 2018-08-29 12:16:24
tags: ['Python', 'Pandas']
---

Recently, I've been using Pandas a lot and have compiled some commonly used functions for future reference.

<!-- more -->


# Data Reading

- Basic reading can be done directly with `read_csv`, which involves several parameters.

```python
import pandas as pd
df = pd.read_csv(file_path, sep="\t", index_col=0, header=0)
# index_col is equivalent to row.names in R; it specifies a column to use as the row labels and removes that column from the data. If not specified, a 0-length numeric label will be generated.
# header functions are different from R's logic; by default, header=0, which means the first row is used as the header. If you don't want a header, use header=None. If you specify another row as the header, rows above that header will be discarded.
```

- When dealing with large datasets to prevent memory overflow, it's necessary to read the file in chunks.

```python
import pandas as pd
# reader is a generator
reader = pd.read_csv(file_path, sep="\t", iterator=True, chunksize=1000)
# iterator / chunksize specifying either will generate an iterator. If you specify chunksize, you don't need to write iterator.

# If reader was specified with chunksize, you can directly call get_chunk()
# To obtain a specific number of rows; otherwise, use get_chunk(num), which I haven't thoroughly tested.
# Note: The row count here refers to the data rows and does not include the header because when generating reader,
# it defaults to using the first row as the header. Therefore, all subsequent generated dfs have headers.
chunk = reader.get_chunk()
df = pd.DataFrame(chunk)
```

- Additionally, Pandas can intelligently recognize whether a file is compressed (by extension) and read/write files with specified delimiters without issues.

{% folding yellow:: Deprecated `read_table` %}

I don't remember from which version onwards, `read_table` was deprecated in favor of `read_csv`, which is said to have better performance for CSV files.

However, it's quite annoying that bioinformatics fields often use `tsv` files or similar formats...

{% endfolding %}

# Data Writing

- My most commonly used formats are `tsv` and `xlsx`.

## Plain Text File

```python
```

## Write to Excel

- The method is very simple:

```python
data.to_excel(FILE, sheet_name='SHEET')
```

- Note that if you want to write multiple data tables to multiple sheets in one go, you can do this:

```python
data.to_excel(FILE, sheet_name=['SHEET1', 'SHEET2'])
```


# Subset Selection

- Similar to R, Pandas supports various subset selection methods. The logic for some selections is similar to R, but there are differences. Generally, whatever can be achieved in Pandas can also be done in R, but the reverse may not be true. Here, I only record methods similar to R's selection logic.

```python
import pandas as pd
df = pd.read_table(file_path, sep="\t")
# Use .loc[row, col] for selecting rows and columns based on labels.
df = df.loc[['a', 'b', 'c'], ['d', 'e', 'f']]
# Note that unlike R, when selecting a single column, there should be no space before the comma. When selecting a single row, you can leave the space after the comma.
df = df.loc[:, ['d', 'e', 'f']]
df = df.loc[['a', 'b', 'c'], ]
# Corresponding to .loc is .iloc, which uses row/column indices instead of labels. Note that it starts from 0.

```

# Row and Column Deletion


# Iteration Processing
