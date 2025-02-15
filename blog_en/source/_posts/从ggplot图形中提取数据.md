---
title: Extracting Outliers from ggplot Boxplots and Removing Them
categories: Script
date: 2020-01-05 17:56:20
tags: ['R', 'ggplot']
---

Today, I encountered a plotting issue and learned a new small trick. I'll record it here.

<!-- Summary -->
<!-- more -->

ggplot provides many commonly used plotting functions, which are quite convenient, and the plots can be drawn in various groupings. However, sometimes I have some special requirements that ggplot cannot meet, so I need to perform some operations manually.

For example, today I need to draw a grouped boxplot, where each group has its own outliers. These outliers will crowd one side of the plot, making it difficult to view, so they need to be removed. Although `geom_boxplot` has parameters for outlier handling, unfortunately, these parameters do not remove the outlier points from the data but simply hide them on the plot. As a result, the crowded plots remain crowded and fail to achieve the desired effect.

The most common solution I found online is to use the `$out` content inside the built-in `boxplot` function to extract outlier values and remove them. Unfortunately, this method cannot provide grouped outlier results; it gives overall outliers instead. If you need to remove outliers by group, you still have to manually split the dataset and process each one individually.

Is there a lazier way?

I found out one myself...

The actual plotting data inside the ggplot object can be extracted:

```r
plot <- ggplot(data=data,mapping=aes(x=data$label, y=data[,tag])) +
    geom_boxplot()
plot_data <- layer_data(plot)
```

The `layer_data` function is one way to extract the data; and `build_plot(plot)$data` has the same function.

`plot_data` stores all the data used by ggplot for plotting. The part we need is stored in `outliers`, and it is arranged in the order shown on the plot. Therefore, as long as we manually specify the factor levels for grouping during the initial plotting, we can obtain the outliers for specific groups from the plot.

Then, utilizing ggplot's automatic handling of `NA` values, we can replace the outlier values in the corresponding column with `NA`, creating a dataset without outliers. Then, we can redraw the plot to solve the problem!

Complete plotting code and plot comparison:

```r
library('ggplot2')
data <- read.csv("plot_data.tsv", sep="\t", stringsAsFactors=F, na.strings = '.')
# Custom levels for grouped data; this must be done subsequently
data$grp <- factor(data$grp, levels=c("Grp1","Grp2"))
# Prepare a vector for mapping
classV <- c(1,2)
names(classV) <- c("Grp1","Grp2")
# First round: use ggplot to get outliers
plot <- ggplot(data=data,mapping=aes(x=data$grp, y=data[,'value'])) +
    geom_boxplot()
print(layer_data(plot))
outlier_data <- layer_data(plot)['outliers']
# Use the obtained values to change outlier values in the original dataset to NA
pdata<-data
for (grp in names(classV)) {
    pdata[pdata$grp==grp,  'value'] <- replace(pdata[pdata$grp==grp, 'value'], pdata[pdata$grp==grp,  'value'] %in% outlier_data[classV[grp], ][[1]], NA)
}
# Prepare a comparison vector; use the changed data to plot
comb_list <- list()
comb <- combn(c("Negative", "Unkown","Weak", "Medium","Strong"),2)
for (id in seq(dim(comb)[2])) {
    comb_list[[id]] <- comb[, id]
}
pplot <- ggplot(data=pdata,mapping=aes(x=pdata$grp, y=pdata[,'value'])) +
geom_boxplot()
```

![display_boxplot.png](https://raw.githubusercontent.com/SilenWang/Gallary/master/display_boxplot.png)