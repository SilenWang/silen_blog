---
title: Registering New Python or R Kernels in Jupyter
categories: Script
date: 2023-01-17 12:53:40
tags: ['jupyter']
---

When testing, it's often necessary to create a new conda environment and then install the jupyter kernel within that environment for use with notebooks. Re-registering new kernels each time can be cumbersome, so here is a record of the process...
<!-- Abstract section -->
<!-- more -->

## Preparing and Registering Python Kernel

- Activate the environment you want to register the kernel in.

```bash
conda activate YOUR_ENV
conda install jupyter ipykernel
python -m ipykernel install --user --name ENV_NAME
```

## Preparing and Registering R Kernel

- Activate the environment you want to register the kernel in.

```bash
conda activate YOUR_ENV
conda install jupyter r-irkernel
```

```r
IRkernel::installspec(name = 'REG_NAME', displayname = 'ENV_NAME')
```
