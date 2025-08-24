---
title: Some Magic Methods in Python
categories: Coding
date: 2019-11-25 23:34:06
tags: ['Python']
---

I previously used `__enter__` and `__exit__`, which I had heard of before. This time, my program development is mainly object-oriented, so I needed to use more methods to achieve some useful features in Python. Here's a small note.

### `__len__`

After defining it, you can use the `len()` function on an instance.

```python
def __len__(self):
    '''
    Define length
    '''
    return len(self.Seq)
```

### `__getitem__`

After defining it, you can perform indexing and slicing operations on an instance using `[]`.

```python
def __getitem__(self, key):
    '''
    Add slice index
    '''
    return self.d[key]
```

### `__str__`

After defining it, when you print an instance using `print()`, it will print the content returned by this method.

```python
def __str__(self):
    return str(self.__dict__)
```

### `__iter__`

After defining it, you can iterate over an object, such as `for i in object`. Note that this method must return an iterator (not just an iterable).

```python
def __iter__(self):
    return (mut for mut in self.muts)
```


### `__eq__` and `__hash__`

These two are together because if you want to automatically remove duplicates when placing instances in a `set()`, both methods need to be defined. Otherwise, the purpose will not be achieved. The `__eq__` method makes objects comparable (`==` and `!=`, `is` cannot be used). Note that the `other` parameter is fixed here; do not change it.

```python
def __eq__(self, other):
    '''
    Ensure automatic removal of duplicates using set characteristics
    '''
    return self.muts == other.muts
```

`__hash__` then essentially makes objects hashable. The hash value used for hashing will be the one provided by this method.

```python
def __hash__(self):
    '''
    Ensure automatic removal of duplicates using set characteristics
    '''
    return hash(self.muts)
```

That's it for now. I'll continue to add more as needed.
```