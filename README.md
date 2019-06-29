
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mipred

The goal of mipred is to calibrate a prediction rule using generalized
linear models or Cox regression modeling, using multiple imputation to
account for missing values in the predictors as described by Mertens,
Banzato and de Wreede (2019) (<https://arxiv.org/abs/1810.05099>).
Imputations are generated using the R package mice without using
outcomes on observations for which the prediction is generated. Two
options are provided to generate predictions. The first is
prediction-averaging of predictions calibrated from single models fitted
on single imputed datasets within a set of multiple imputations. The
second is application of the Rubin’s rules pooled model. For both
implementations, unobserved values in the predictor data of new
observations for which the predictions are derived are automatically
imputed. The present version of the package is preliminary (development)
and has only been checked to support binary-outcome logistic regression
for now. We are working to expand the functionality to non-binary and
survival
outcomes.

## Installation

<!-- You can install the released version of mipred from [CRAN](https://CRAN.R-project.org) with: -->

You can install the current version into R from GitHub using devtools:

``` r
devtools::install.github("BartJAMertens/mipred")
```

You may need to install and load the devtools package first before using
the above command. See the book “R packages” (online version) by Hadley
Wickham, chapter “Git and Github”.

## Example

Please refer to the example included with the package. The package also
includes a vignette which documents use for binary outcome
data.

<!-- This is a basic example which shows you how to solve a common problem: -->

<!-- ```{r example} -->

<!-- ## basic example code -->

<!-- ``` -->

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->

<!-- ```{r cars} -->

<!-- summary(cars) -->

<!-- ``` -->

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. -->

<!-- You can also embed plots, for example: -->

<!-- ```{r pressure, echo = FALSE} -->

<!-- plot(pressure) -->

<!-- ``` -->

<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub! -->
