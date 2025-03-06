---
title: R Training Notes
date: 2018-07-29 00:35:55
tags: ['R']
categories: Script
---

Graduate training on R at the Guangdong Provincial CDC was my first serious programming course.

<!-- more -->

# R Training Notes Summary

## Overview
***
The training lasted three days and covered five main topics:

- Introduction to R language and basic concepts, operations
- Using R to draw statistical graphs (graph & ggplot2)
- Introduction to time series analysis and the use of Generalized Linear Models (GLM) and Generalized Additive Models (GAM)
- Use of Distributed Lag Nonlinear Models
- Meta-analysis in R

## Introduction to R Language and Basic Concepts, Operations
***

### R Features

- Pros: Free, rich features, fast updates, process visualization, beautiful plots
- Cons: Difficult to master, inconsistent, requires time to learn, slow (can be optimized with specific packages)

### Software Installation

- R and RStudio can be downloaded from the official website
- Packages can be installed manually or through the software itself
- Functions in R need to load the corresponding package before use, and help can only be viewed after loading.

### Software Interface Introduction (omitted)

### Common R Syntax

- Case-sensitive
- Statements are composed of functions and assignments
- All symbols must be half-width English
- `#` is used for comments, `#text####` can be indexed (script mode)
- No special symbol at the end of each line; if a statement is incomplete, the program will continue to execute on the next line until it is complete
- `<-` is the assignment operator, e.g., `x<-3`, `=` can replace it
- `==` is the logical equality operator, `&` for AND, `|` for OR, `>=` greater than or equal to, `<=` less than or equal to
- Single and double quotes are mostly interchangeable
- There are two equivalent ways to express paths: `C:\\R\\Data` is equivalent to `C:/R/Data`; R has poor support for Chinese,尽量 use full English paths
- `$` can be used to extract variables within an object; after extraction, the variable can be assigned or modified.
- `[]` can be used to extract subsets of a data frame
- `!` is the negation symbol

### Setting Working Directory & Data Reading/Writing

- Directory:
`getwd("path")/setwd("path")`
- Data Import:
    * Excel (only supports .csv):
    `read.csv("path+filename, header= , sep=",")`
    * txt:
    `read.table("path+filename", header= , )`
    * SPSS (.sav):
    `read.spss("path+filename", use.value.labels= )`
    * SAS (.sas7bdat):
    `read.sas7bdat("path+filename")`
    * R (.rda):
    `load("path+filename")`
    * Direct read:
    `read.table(textConnection("content"), header=, sep=)`
- Data Export
    * R:
    `save(对象, file="path+filename")`
    * Excel (.csv):
    `write.csv(对象, file="path+filename")`
    * txt:
    `write.table(对象, file="test.txt",sep = "delimiter")`

### Viewing Objects (Data)

- Multiple functions can be used to view already loaded (created) objects, each function has its own focus.
    * `View(对象)` opens a window on the left to display object/data information;
    * `head(对象, 条目数)` displays the first n rows of data in the console; if no number is set, it defaults to 6 rows;
    * `tail(对象, 条目数)` displays the last n rows of data in the console; if no number is set, it defaults to 6 rows;
    * `names(对象)` displays all variable names (headers) of the object in the console.
    * `str(对象)` displays the number of entries and variables of the object, as well as each variable's name, type, and the first 10 values.
    * `dim(对象)` displays the number of rows and columns of the object

### Statistical Description of Objects (Data)

Function basic form|Function purpose|Notes
-------------|-------------|------
`summary(对象)`|Displays maximum, minimum, arithmetic mean, median, quartile range for numerical variables; <br>Counts for categorical variables; <br>If missing values exist, displays the number of missing values.|Use `rm.na=` parameter to automatically remove missing values when calculating corresponding statistics.

### Simple Operations on Objects (Data)

Function expression|Operation effect|Notes
----------------|-------------|-------
`names(data)<-c("1","2","3")`|Renames the first three variables of data frame data to 1, 2, 3|If there are more than three variables in the data frame, the names of the variables after the third will be lost; if too many names are assigned, it cannot be executed
`data=rename(data,c("old1"="new1", "old2"="new2"))`|Replaces specified variable names in data frame data|rename() is a function from the plyr package
`as.date(data)`<br>`as.numeric(data)`<br>`as.character(data)`<br>`as.factor(data)`<br>etc.|Converts specified variables to types, after conversion, can be assigned to an object|
`subset(data, select=c(v1,v2))`<br>`subset(data, v1==n)`<br>`data[row,colum]`|Extracts specified rows/columns from data frame|Whether the two methods are completely equivalent is unknown; `t[]` defaults to taking rows
`data$v=data$v1+test$v2`|Performs calculations on existing data in data frame and writes it as a new variable|
`substr(strData,1,3)`|Extracts a specified number of characters from a string variable|
`mutate(data,)`|Simplifies code by performing multiple operations simultaneously on data frame data|The main function of mutate() is to simplify code
`join(data1, data2, by= ,type=)`|Merges different data frames according to certain standards|join() function does not require sorting of the table

## Using R to Draw Statistical Graphs

***
- Two common graph drawing packages in R are:
  * graph (built-in): Basic plotting
  * ggplot2 (needs installation): Advanced plotting, high customization, beautiful finished products
