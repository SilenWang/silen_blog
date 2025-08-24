---
title: Python Command Line Parsing
categories: Coding
date: 2018-10-22 12:12:50
tags: ['Python', 'argparse', 'Argument Parsing']
---

I have been using `argparse` for command line parsing, but it wasn't until recently that I discovered it can be used to create subcommands.

<!-- more -->

## General Argument Parsing

General argument parsing is quite simple. You just instantiate an `argparse`, and then call the object's add argument method.

```python
import argparse
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(
    prog='program',
    description='Some description',
    formatter_class=RawTextHelpFormatter
)
parser.add_argument('-s',
                    help="description for args",
                    required=False,
                    type=str)
```

## Subcommand Creation

Subcommands need to first create a subparser object based on the main argument parser (above is `parser`), and then call this object's method to add subcommands.

```python
sub_parser = parser.add_subparsers(title='sub-commands',
                                    description='Some description',
                                    help='Some help description')
sub_parser_a = sub_parser.add_parser('a', help='Some help description')
sub_parser_a.add_argument('-s', 
                        help="description for args",
                        type=str)
```

This way, you can define subprograms under the main program. `argparse` also provides a method to specify a function to run when a specific subcommand is called. This function will receive all parameter values under the subcommand as a namespace and can start executing the subprogram after parsing the parameters.

```python
sub_parser_a.set_defaults(func=subFun)

def countFun(args):
    args = vars(args)
    print(args)
```

## Other

`argparse` can also implement multiple subprograms and inherit main program parameters for subprograms. This has not been used yet, so I will update it later.