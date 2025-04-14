---
title: Multiprocessing in Python
categories: Script
date: 2019-02-12 23:45:59
tags: ["Python", "Multiprocessing", "Parallel Computing"]
---

R provides convenient multiprocessing capabilities, and Python has similar functionality.
<!-- Excerpt -->
<!-- more -->
First, let's quote from [Liao Xuefeng's Python Tutorial](https://www.liaoxuefeng.com/wiki/0014316089557264a6b348958f449949df42a6d3a2e542c000/0014319272686365ec7ceaeca33428c914edf8f70cca383000) to explain processes and threads:

> For an operating system, a task is a Process. For example, opening a browser starts a browser process, opening Notepad starts a Notepad process, opening two Notepads starts two Notepad processes, and opening Word starts a Word process.
>
> Some processes do more than one thing at the same time. For example, Word can perform typing, spell checking, and printing simultaneously. Within a process, to perform multiple "subtasks" at the same time, we call these "subtasks" Threads.
>
> Since each process must do at least one thing, a process has at least one thread. Of course, complex processes like Word can have multiple threads that execute simultaneously. The execution of multiple threads is similar to multiple processes - the operating system quickly switches between threads, making each thread run alternately for short periods, appearing to execute simultaneously. True simultaneous execution of multiple threads requires a multi-core CPU.
>
> All the Python programs we've written before are single-task processes with only one thread. How do we execute multiple tasks simultaneously?
> There are two solutions:
> One is to start multiple processes, where each process has only one thread, but multiple processes can execute multiple tasks together.
>
> Another method is to start one process and multiple threads within that process, so multiple threads can execute multiple tasks together.
>
> There's also a third method - start multiple processes, each with multiple threads - but this model is more complex and rarely used in practice.

The tasks I need to execute don't involve long IO waits, so I didn't choose Python's multithreading or coroutines, but instead tried multiprocessing.

The main code is as follows:

```python
from multiprocessing import Pool

def func_to_process(args):
    do someting
    return result


def main():
    sub_procs = Pool(10)  # 10 is the number of child processes to use
    results = [] # To store results
    for obj in iterable:
        result = sub_procs.apply_async(func_to_process, (args))
        # Submit task to process pool, apply_async means submit and continue
        # Will block when pool is full
        results.append(result)
    sub_procs.close() # Close pool, no new tasks can be submitted
    sub_procs.join() # Wait for all processes to complete
    results = [result.get() for result in results] # Use get() to retrieve all results


if __name__ == "__main__":
    main()
```

There are two issues with this code:

1. This approach splits the file by lines with each line as a process, which is too fragmented. Especially when the task has lengthy preparation steps, this actually reduces efficiency (single process prepares once, while multiprocess prepares many times)
2. The results may be out of order. If the task requires maintaining order, the results need to be sorted again

For issue 1, after some research, it might be possible to use some functions from `itertools` to implement quantitative splitting of the original file. Will update tomorrow.

For issue 2, we could use the corresponding map method instead of apply. This requires further testing.
