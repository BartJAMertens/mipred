#' Cross-validation of generalized linear model prediction using multiple imputation - prediction-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param nimp Number of imputations for each observation
#' @param folds Number of folds defined in newdata
#' @param miop Mice options

.glm.mipred.cmb1.cv <-
  function(formula, family, dataset, nimp, folds, miop) {

    n <- nrow(dataset)

    Xb  <-
      Pred <-
      matrix(NA, nrow = n, ncol = nimp)   # Initialize predictor matrices for new observations

    for (m in 1:nimp) {
      suppressWarnings(folddef <-
          split(sample(n, n, replace = F), as.factor(1:folds)))

      for (k in 1:folds) {
        folddef[[k]] <-
          sort(folddef[[k]])    # Sort the ids of the k-th fold

        datanew <-
          dataset             # Make a data copy of the original data to keep the latter unaffected
        datanew[folddef[[k]], 1] <-
          NA      # Remove the response from the validation data partition

        # Compute MI (once) in the whole dataset, run an imputation with m = 1, use previously generated options

        imp_data <-
          impute(datanew, miop, 1, miop[["seed"]][k + (m - 1) * folds])

        data_compl <-
          complete(imp_data, 1)  # Select the completed data

        # Fit a model on the completed calibration data only, remove the validation data partition
        fit <-
          glm(formula, family = family, data = data_compl[-folddef[[k]], ])

        coefs <- fit$coefficients # Save the model coefficients
        X_m <-
          model.matrix(formula,
            data = data_compl[folddef[[k]],],
            drop.unused.levels = FALSE) # Define the regression matrix in validation partition, prevent dropping of unused factor levels

        Xb[folddef[[k]], m] <-
          X_m[, names(fit$coefficients)] %*% coefs # Compute X*b only for the validation sample

        Pred[folddef[[k]], m] <-
          predict.glm(fit, data_compl[folddef[[k]], ], type = "response") # This will generate error if more factor levels than calibration set
      }
      print(m)
    }
    output <- list(pred = Pred, linpred = Xb)
    output
  }

#' Cross-validation of generalized linear model prediction using multiple imputation - Rubin's rule coefficient-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param folds Number of folds defined in data
#' @param nimp Number of imputations for each observation
#' @param miop Mice options

.glm.mipred.cmb2.cv <-
  function(formula, family, dataset, nimp, folds, miop) {

    attributes(dataset)$terms <- NULL
    n <- nrow(dataset)

    modeldim <-
      dim(model.matrix(formula , data = dataset))[2] # Length of model coefficient vector

    Xb <-
      Pred <-
      matrix(NA, nrow = n, ncol = nimp) # Initialize predictor matrices for new observations

    coefs <-
      matrix(NA, nrow = folds, ncol = modeldim) # Only if we want to save the results - Rubin's rule pooled coefs
    suppressWarnings(folddef <-
        split(sample(n, n, replace = F), as.factor(1:folds))) # folddef is constant for Rubin's rule application

    for (k in 1:folds) {
      folddef[[k]] <- sort(folddef[[k]])    # Sort the id of the k-th fold
      datanew <-
        dataset   # Make a copy of the original data to keep the latter unaffected

      datanew[folddef[[k]], 1] <-
        NA # Remove the response from the validation partition
      X_m <-
        array(NA, dim = c(length(folddef[[k]]), modeldim, nimp)) # Initialize array to store imputed validation folds across imputations

      # Compute MI with m=M in the whole dataset, run an imputation with m = M, use previously generated options
      imp_data <- impute(datanew, miop, nimp, miop[["seed"]][k])
      coef_m <-
        matrix(NA, nrow = nimp, ncol = modeldim) # Initilize coefs matrix prior to Rubin's rule averaging

      for (m in 1:nimp) {
        data_compl <- complete(imp_data, m)  # Select the completed data

        # Fit a Logistic Cox model on the completed training data, remove the validation data portion
        fit <-
          glm(formula , family = family, data = data_compl[-folddef[[k]], ])
        X_m[, , m]  <-
          model.matrix(formula,
            data = data_compl[folddef[[k]],],
            drop.unused.levels = FALSE) # Save the mth completed validation model matrix
        coef_m[m, ] <-
          fit$coefficients       # Save the coefficients across all imputations
      }
      coefs[k, ] <-
        apply(coef_m, 2, mean)     # Pool the coefficients and save them
      Xb[folddef[[k]], ] <-
        apply((X_m) * ((rep(
          1, length(folddef[[k]])
        )) %*% t(coefs[k,])) %o% (rep(1, nimp)), c(1, 3), sum) # Apply pooled coefficients to all imputed validation data matrices and save linear predictors

      fit$coefficients <-
        coefs[k,] # replace coefficient vector with pooled coefficients

      if (!any(is.na(dataset[folddef[[k]],]))) {
        # complete record - so all predictions will be the same within individual
        Pred[folddef[[k]],] <-
          matrix(predict.glm(fit, data_compl[folddef[[k]], ], type = "response"), ncol =
              1) %*% rep(1, nimp)
      } else {
        for (m in 1:nimp) {
          data_compl <- complete(imp_data, m)
          Pred[folddef[[k]], m] <-
            predict.glm(fit, data_compl[folddef[[k]], ], type = "response") # This will generate error if more factor levels than for calibration sets
        }
      }
      print(k)
    }
    # Pred2 <-
    #   .expit(Xb) # Compute the probability for each individual in the validation set and save the predictions
    # if (sum(abs(Pred - Pred2)) > 1e-6)
    #   browser()
    output <- list(pred = Pred, linpred = Xb)
    output
  }
