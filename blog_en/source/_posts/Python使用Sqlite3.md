---
title: Python使用Sqlite3
categories: Script
date: 2018-08-26 23:54:24
tags: ['Python', 'sqlite']
---

To prepare for future operations using Sqlite3, I recently tried using it to operate on databases at work. I decided to use sqlite3 because it seemed convenient. However, I found that not all SQL statements are supported by sqlite3, and its free-form shell commands cannot be called directly by the Python sqlite module...

<!-- more -->

## Simple Usage:

```python
import sqlite3   # Import the module
conn = sqlite3.connect('test.db')   # Connect to the database (it will be created if it doesn't exist)
cursor = conn.cursor()    # Create a pointer from the database object
cursor.execute('sql cmd')   # Pass SQL statements to the database, note that operations are not immediately saved/synchronized to the database
cursor.close()   # Close the pointer
conn.commit()   # Confirm changes
conn.close()   # Close the connection
```

- Usage Notes:
    * Only standard SQL statements are supported. Since sqlite is being used, unsupported SQL syntax cannot be used, and sqlite-specific shell commands cannot be used either because these commands are essentially executed in the shell, not as SQL statements.
    * If you must use shell commands, you can use the os module or another module that can call the shell, effectively executing shell commands.
    * Not all operations need to be completed through the pointer object. Some can be done using the database object itself, which can save a few lines of code. See [the introduction on the Coder's Tutorial website](http://www.runoob.com/sqlite/sqlite-python.html) for more details.
    * Remember to confirm changes (`obj.commit()`), otherwise the content in the database will not change.
    * I tried using the with statement to connect, but it doesn't automatically close the connection outside the with statement scope (`obj.close()`), so you still need to write it explicitly.

## Calling Shell Commands

```python
import os, textwrap

# This is where textwrap is used to adjust shell command indentation for better readability when outputting. It's not necessary.
def cmd_gen(cmd_str):
    cmd_str = textwrap.dedent(cmd_str)
    cmd_str = cmd_str.strip()
    return cmd_str


db_file = 'test.db'
index_file = 'index.txt'

# Generate the command, this is the only valid method I found online. The specific meaning will be explained later...
cmd = cmd_gen('''
    sqlite3 {db_file} << EOF
    .separator "\\t"
    .import {index_file} idx_tab
    EOF
    '''.format(db_file=db_file, index_file=index_file))
os.system(cmd)
```
```