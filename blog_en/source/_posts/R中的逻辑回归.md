---
title: Logistic Regression in R
categories: Others
date: 2019-03-16 14:07:19
tags: ["R"]
---

I didn't expect to do regression again after four years... this time using R instead of SPSS. Previously, I was doing statistical analysis, but now... it's machine learning.
<!-- Abstract part -->
<!-- more -->

## Introduction
Logistic (logit) regression is a type of generalized linear model that has the characteristic of converting a categorical dependent variable into a continuous variable by comparing all categories to one reference category. In logistic regression, the coefficients of the variables are taken as `exp(x)` after which they represent the OR value for this variable, meaning that for every unit increase in this variable, the probability of being classified as 1 increases by `exp(x)`.

## Model Construction
This time, my dependent variable is binary, so I can use the basic package's `glm()` function. The specific model construction uses:

- `family=binomial(link='logit')` to specify calling a binary logistic regression model

```r
model <- glm(y ~ x1 + x2, family=binomial(link='logit'), data=fit_data)
summary(model) # Display model information
```

When performing regression, it's difficult to guarantee that the first regression will yield the best results. Among the included covariates, there may be variables with coefficients too small (meaning they have little effect), or statistical tests are not significant (which may not have an effect in a resampling). Therefore, we need to change variables multiple times for regression.

In R, `step()` can help automatically perform regression and remove variables that contribute little and are not significant:

```r
model.step <- step(model)
summary(model.step) 
```

Of course, this is just an automated reference. Ultimately, how the model is determined still depends on specific judgment by a person based on actual circumstances.

## Model Prediction Situation Estimation

After building the model, you can use `predict()` to predict test datasets and evaluate the model. Then, use some third-party packages through ROC curves to determine fitting situations and select final classification thresholds.

```r
test_data$prob <- predict(model.step, test_data, type = "response") # response means getting a probability between 0~1, not a constant

# Use the pROC package
library(pROC)

obj_roc <- roc(test_data$Real, test_data$prob) # Real is the correct answer, prob is the predicted probability
plot(obj_roc, print.auc=TRUE, auc.polygon=TRUE, grid=c(0.1, 0.2), # Make a plot
     grid.col=c("green", "red"), max.auc.polygon=TRUE,
     auc.polygon.col="skyblue", print.thres=TRUE)
# The pROC package gives an ROC curve graph with sensitivity (true positive) on the y-axis and specificity (true negative) on the x-axis, not the more common false positive rate, so the x-axis is 1->0 instead of 0-1
```

## Small Tips

The model object itself contains detailed model information, but directly taking information from the model does not match what `summary()` sees. For example, the coefficient table, `summary()` gives a table that includes coefficient values, standard errors, and statistical test P-values, but directly taking it from the model will only have coefficient values.

To take `summary()` information, you need to take it from the generated object or assign it first:

```r
# First assignment
model_info <- summary(model)
model_info$coefficients
# Directly take
summary(model)$coefficients
```

## Blog Article Learning

I found a few blog posts on using R for logistic regression to supplement model evaluation methods.

### [Logistic Regression](http://r-statistics.co/Logistic-Regression-With-R.html)

- Balance of binary variables: In logistic regression, if the two categories of the dependent variable are balanced, the classification results will be ideal. Therefore, before performing regression, you should check the situation of the two categories: `table(your_data$val)`. If they are not balanced, theoretically, you should sample the larger category to match the number of the smaller category.
- Multicollinearity detection: Logistic regression and linear regression require routine detection of multicollinearity issues. The blog recommends `VIF<4`.

### [How to perform a Logistic Regression in R](https://www.r-bloggers.com/how-to-perform-a-logistic-regression-in-r/)

- $R^2$ calculation: Logistic regression cannot calculate $R^2$. The article suggests calculating the *McFadden $R^2$ index* as an alternative. More information about this metric can be found at [R squared in logistic regression](http://thestatsgeek.com/2014/02/08/r-squared-in-logistic-regression/).

### [Logistic regression on biased data](https://datascience.stackexchange.com/questions/12234/logistic-regression-on-biased-data)

- Evaluation metrics selection in imbalanced classification: When the categories are imbalanced, `accuracy` is not a good evaluation metric. It is recommended to use `F1-score`.
- Solutions for imbalanced classification include two types: [oversampling the minority class](https://www.marcoaltini.com/blog/dealing-with-imbalanced-data-undersampling-oversampling-and-proper-cross-validation) and fixing the model by altering the [hyperplane (SVM)]() or [changing priors (Bayes)]().
- According to the article, oversampling the minority class and undersampling the majority class can improve detection results to some extent, but they are not good enough for testing data.

```r
library(pscl)
pR2(model)
```

### Other Materials

- [Logistic Regression Algorithm](https://blog.mythsman.com/2016/01/28/1/)
- [Cross-validation idea in machine learning](https://blog.mythsman.com/2016/02/02/1/)
- [Dealing with unbalanced data in machine learning](https://www.r-bloggers.com/dealing-with-unbalanced-data-in-machine-learning/)
```