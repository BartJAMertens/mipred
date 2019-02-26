# mipred
mipred
======

The goal of mipred is to calibrate a prediction rule using generalized linear models or Cox regression modeling, using multiple imputation to account for missing values in the predictors as described by Mertens, Banzato and de Wreede (2019) (<https://arxiv.org/abs/1810.05099>). Imputations are generated using the package mice without using outcomes on observations for which the prediction is generated. Two options are provided to generate predictions. The first is prediction-averaging of predictions calibrated from single models fitted on single imputed datasets within a set of multiple imputations. The second is application of the Rubin's rules pooled model. For both implementations, unobserved values in the predictor data of new observations for which the predictions are derived are automatically imputed. The present version of the package is preliminary (development) and has been checked to only support binary-outcome logistic regression. We are working to expand the functionality to non-binary and survival outcomes.

Installation
------------

<!-- You can install the released version of mipred from [CRAN](https://CRAN.R-project.org) with: -->
You can install the current version from GitHub using devtools:

``` r
devtools::install.github("BartJAMertens/mipred")
```

Example
-------

Please refer to the example included with the package
