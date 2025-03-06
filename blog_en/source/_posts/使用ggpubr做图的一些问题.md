---
title: R Plotting Issues Record
categories: Script
date: 2022-03-13 03:12:07
tags: ['R', 'ggplot', 'ggpubr', 'survminer']
---

Recently, while helping a colleague create plots, I uncovered some new issues and solutions related to `ggpubr` and survival analysis packages. Here are the records.

<!-- Abstract section -->
<!-- more -->

## survminer Plot Saving Issue

The `survminer` package includes a convenient function called `ggsurvplot`, which allows for easy creation of survival curves with group significance tests, as well as frequency tables below the plot. However, this function returns a `ggsurvplot` object rather than a `ggplot()` object. This object contains `survObj$plot` for the survival curve and `survObj$table` for the frequency table. Therefore, if you need to adjust or modify the plot using `ggplot` syntax, you should operate on the objects within these two components.

## ggpubr Still Has Bug with Significance Test Annotations in Faceted Plots

As of December last year, there is still a bug in `ggpubr` where box plots (other facets may also have this issue) append significance test annotations incorrectly after faceting. The data and annotation positions are incorrect. You need to manually draw individual plots and then combine them.

## Table Annotations Inside the Plot

To be filled in later.

## Multiple Ways to Implement Significance Annotations

Due to a complex requirement from a client, where significance tests should be based on a reference group with two groups instead of one, `ggpubr`'s built-in functions cannot meet this need. I spent considerable time researching how to manually implement it, and different methods for various requirements are recorded here:

To be filled in later.
