---
title: Using the parallel package in R
categories: Script
date: 2019-01-01 02:04:56
tags: ['R', 'Parallel Computing']
---

Recently, I wrote a script to draw many images using ggplot. Although all the places where loops are needed have already been replaced with `apply`, it still can't keep up when drawing hundreds of images at once. So, I fiddled around with `parallel` and managed to parallelize the plotting part to speed things up.
<!-- Abstract section -->
<!-- more -->
Before using Python, I was familiar with multithreading and multiprocessing. However, in R... it doesn't seem to emphasize threads and processes as much; instead, it focuses on parallel computing. The only relevant resources I found in Chinese were `parallel` and `foreach`. Of course, I used `parallel`, but I haven't studied the differences between the two...

Using it is actually quite simple. First, you register a cluster, then you throw your tasks at the cluster to run (in function form), and finally, you stop the cluster. The "cluster" here doesn't refer to a cluster server; it's just what English calls 'cluster', which I don't know how to translate well... I'll just call it that for now. Here is the specific implementation code:

```R
# Prepare the plotting function
plot_process <- function (data_list, args) {
    library(ggplot2)
    plot <- ggplot() + geom_line()
    ggsave("plot.png", plot)
}

cl <- makeCluster(core_num) # Register
parLapply(cl, data_list, plot_process, args) # The first parameter specifies the cluster, and the rest are the same as lapply
stopCluster(cl) # Stop
# Since I'm only processing data separately and plotting, there's no need for the merging step that many tutorials have.
```

So it's actually quite simple... The only slightly troublesome issue is that for each cluster node, you need to load the variables, packages, and functions separately. For variables, there are specific functions to pass them; for packages, you can write the loading statements inside `plot_process`. However, for custom functions, there isn't a direct way to do it...
After searching for a long time without finding a solution, I suddenly had an idea of passing the functions as arguments to `plot_process`, and it worked...

After writing this script, I also saw another solution: write these custom functions in another file and then import them using `source()` inside `plot_process`. This should be a better approach, so I'll try it out later....
