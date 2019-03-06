#' Cross-validation of binary prediction using multiple imputation - prediction-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param folds The fold definition
#' @param nimp Number of imputations for each observation

.glm_mipred_cmb1_cv <- function(formula, family, dataset, nimp, folds){

  formulapredict <- formula
  formulapredict[[2]] <- NULL # Use delete.response?
  n <- nrow(dataset)

  Xb  <- matrix(NA, nrow = n, ncol = nimp)   # Linear predictor for new observations

  for (m in 1:nimp){

    folddef <- split(sample(n, n, replace=F), as.factor(1:folds))

    for (k in 1:folds){

      folddef[[k]] <- sort(folddef[[k]])    # Sort the ids of the k-th fold

      datanew <- dataset             # Make a data copy of the original data to keep the latter unaffected
      datanew[folddef[[k]], 1] <- NA      # Remove the response (Status) from the validation data partition

      # Compute MI (once) in the whole dataset, run an imputation with m = 1
      imp_data <- mice(datanew, m = 1, maxit = 50, printFlag = FALSE)
      data_compl <- complete(imp_data, 1)  # Select the completed data

      # Fit a Logistic Cox model on the completed calibration data only, remove the validation data partition completely
      fit <- glm(formula, family = family, data = data_compl[-folddef[[k]], ]) # Default family=binomial("logit")
      coef <- fit$coefficients # Save the model coefficients
      X_m <- model.matrix(formulapredict, data = data_compl[folddef[[k]],]) # Define the regression matrix in validation partition
      Xb[folddef[[k]], m] <- X_m %*% coef # Compute X*b only for the validation sample
    }
    print(m)
  }
  pred <- .expit(Xb)
  output <- list(pred = pred, linpred = Xb)
  output
}

#' Cross-validation of binary prediction using multiple imputation - Rubin's rule coefficient-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param folds The fold definition
#' @param nimp Number of imputations for each observation

.glm_mipred_cmb2_cv <- function(formula, family, dataset, nimp, folds){

  formulapredict <- formula
  formulapredict[[2]] <- NULL   # Use delete.response?
  n <- nrow(dataset)

  modeldim <- dim(model.matrix(formula , data = dataset))[2] # Length of model coefficient vector

  Xb <- matrix(NA, nrow = n, ncol = nimp)   # Pred <- matrix(NA, nrow=n, ncol=nimp)
  coefs <- matrix(NA, nrow = folds, ncol = modeldim) # Only if we want to save the results - Rubin's rule pooled coefs
  folddef <- split(sample(n, n, replace = F), as.factor(1:folds)) # folddef is constant for Rubin's rule application

  for (k in 1:folds){

    folddef[[k]] <- sort(folddef[[k]])    # Sort the id of the k-th fold
    datanew <- dataset   # Make a copy of the original data to keep the latter unaffected

    datanew[folddef[[k]], 1] <- NA # Remove the response from the validation partition
    X_m <- array(NA, dim = c(length(folddef[[k]]), modeldim, nimp)) # Initialize array to store imputed validation folds across imputations

    # Compute MI with m=M in the whole dataset, run an imputation with m = M
    imp_data <- mice(datanew, m = nimp, maxit = 50, printFlag = FALSE)

    coef_m <- matrix(NA, nrow = nimp, ncol = modeldim) # Initilize coefs matrix prior to Rubin's rule averaging

    for (m in 1:nimp){

      data_compl <- complete(imp_data, m)  # Select the completed data

      # Fit a Logistic Cox model on the completed training data, remove the validation data portion
      fit <- glm(formula , family=family, data = data_compl[-folddef[[k]], ])
      X_m[, , m]  <- model.matrix(formulapredict, data = data_compl[folddef[[k]],]) # Save the mth completed validation model matrix
      coef_m[m, ] <- fit$coefficients       # Save the coefficients across all imputations
    }
    coefs[k, ] <- apply(coef_m, 2, mean)     # Pool the coefficients and save them
    Xb[folddef[[k]], ] <- apply((X_m)*((rep(1,length(folddef[[k]])))%*%t(coefs[k,]))%o%(rep(1,nimp)),c(1,3),sum) # Apply pooled coefficients to all imputed validation data matrices and save linear predictors
  }
  pred <- .expit(Xb) # Compute the probability for each individual in the validation set and save the predictions

  output <- list(pred = pred, linpred = Xb)
  output
}







