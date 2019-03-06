#' Expit function converting odds to probability
#'
#' @param x Probability vector

.expit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Binary prediction using multiple imputation - prediction-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param newdata A dataframe containing observations to be predicted
#' @param nimp Number of imputations for each observation
#' @param folds Number of folds defined in newdata
#' @param miop Mice options

.glm_mipred_cmb1 <-
  function(formula,
    family,
    dataset,
    newdata,
    nimp,
    folds,
    miop) {
    formulapredict <- formula
    formulapredict[[2]] <- NULL # use delete.response?
    newn <- nrow(newdata)

    newdata <-
      base::cbind(NA, newdata) # append column vector of NA for absent response in validation set
    names(newdata)[[1]] <-
      names(dataset)[[1]] # rename to allow rbind below

    # Define some matrices we need later
    Xb <-
      matrix(NA, nrow = newn, ncol = nimp)   # linear predictor for new observations
    Pred <-
      matrix(NA, nrow = newn, ncol = nimp)  # predicted probabilities for new observations
    mioppredictorMatrix <- is.null(miop$predictorMatrix)
    miopblocks <- is.null(miop$blocks)
    miopformulas <- is.null(miop$formulas)

    # check predictorMatrix, blocks and formulas
    if (mioppredictorMatrix & miopblocks & miopformulas) {#browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (!mioppredictorMatrix & miopblocks & miopformulas) {#browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, predictorMatrix=miop$predictorMatrix)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (mioppredictorMatrix & !miopblocks & miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, blocks=miop$blocks)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (mioppredictorMatrix & miopblocks & !miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, formulas=miop$formulas)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (!mioppredictorMatrix & !miopblocks & miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, predictorMatrix=miop$predictorMatrix, blocks=miop$blocks)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
  if (!mioppredictorMatrix & miopblocks & !miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, predictorMatrix=miop$predictorMatrix, formulas=miop$formulas)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (mioppredictorMatrix & !miopblocks & !miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, blocks=miop$blocks, formulas=miop$formulas)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }
    if (!mioppredictorMatrix & !miopblocks & !miopformulas) {browser()
      # Run an imputation with zero iterations to create the predictor matrix
      imp_init <- mice(dataset, m = 1, maxit = 0, predictorMatrix=miop$predictorMatrix, blocks=miop$blocks, formulas=miop$formulas)
      miop$predictorMatrix <- imp_init$predictorMatrix
      miop$blocks <- imp_init$blocks
      miop$formulas <- imp_init$formulas
    }

#browser()
    for (m in 1:nimp) {
      folddef <- split(sample(newn, newn, replace = F), as.factor(1:folds))

      for (k in 1:folds) {
        folddef[[k]] <-
          sort(folddef[[k]])    # Sort the ids of the k-th fold
        combdat <-
          base::rbind(newdata[folddef[[k]],], dataset)  # case-by-case predictions for new observations
        lengthfold<-length(folddef[[k]])
#browser()
        # Compute MI (once) in the combined dataset, run an imputation with m = 1
        imp_data <-
          mice(
            dataset,
            m = 1,
            method = miop[["method"]],
            #predictorMatrix = miop[["predictorMatrix"]],
            where = miop[["where"]],
            blocks=miop[["blocks"]],
            visitSequence =  miop[["visitSequence"]],
            formulas=miop[["formulas"]],
            blots = miop[["blots"]],
            post = miop[["post"]],
            defaultMethod = miop[["defaultMethod"]],
            maxit = miop[["maxit"]],
            printFlag = miop[["printFlag"]],
            seed = miop[["seed"]][k+(m-1)*folds],
            data.init = miop[["data.init"]]
          )
        #print(k+(m-1)*folds)
        #print(miop[["seed"]][k+(m-1)*folds])
        data_compl <-
          complete(imp_data, 1)  # select the completed model data

        # Fit a Logistic Cox model on the completed training data, remove the validation data portion
        fit  <-
          glm(formula, family = family, data = data_compl[-(1:lengthfold),]) # default family=binomial("logit")
        coef <- fit$coefficients # Save the model coefficients

        X_m <-
          model.matrix(formulapredict, data = data_compl[(1:lengthfold), ])   # Define the regression matrix and remove status next
        Xb[folddef[[k]], m] <-
          X_m %*% coef            # Compute X*b only for the validation sample

        Pred[folddef[[k]], m] <- predict.glm(fit,data_compl[(1:lengthfold), ],type="response")
      }
    #print(m)
    }

    # Compute the probability for each individual in the validation set
    #Pred <- .expit(Xb)
#browser()
    output <- list(pred = Pred, linpred = Xb)
    output
  }

#' Binary prediction using multiple imputation - Rubin's rule coefficient-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param newdata A dataframe containing observations to be predicted
#' @param nimp Number of imputations for each observation
#' @param folds Number of folds defined in newdata
#' @param miop Mice options

.glm_mipred_cmb2 <-
  function(formula,
    family,
    dataset,
    newdata,
    nimp,
    folds,
    miop) {
    formulapredict <- formula
    formulapredict[[2]] <- NULL   # use delete.response?
    newn <- nrow(newdata)

    modeldim <-
      dim(model.matrix(formula , data = dataset))[2] # length of model coefficient vector
    newdata <-
      base::cbind(NA, newdata) # append column vector of NA for absent response in validation set
    names(newdata)[[1]] <-
      names(dataset)[[1]] # rename to allow rbind below

    # Define some matrices we need later
    coef_m <- matrix(NA, nrow = nimp, ncol = modeldim)
    X_m <- matrix(NA, nrow = nimp, ncol = modeldim)
    Xb <- matrix(NA, nrow = newn, ncol = nimp)
    Pred <- matrix(NA, nrow = newn, ncol = nimp)

    coefs  <-
      matrix(NA, nrow = folds, ncol = modeldim) # only if we want to save the results
    folddef <-
      split(sample(newn, newn, replace = F), as.factor(1:folds))

    for (k in 1:folds) {
      folddef[[k]] <-
        sort(folddef[[k]])    # Sort the ids of the k-th fold
      lengthfold<-length(folddef[[k]])
      combdat <-
        base::rbind(newdata[folddef[[k]],], dataset)  # case-by-case predictions for new observations
      X_m <- array(NA, dim = c(lengthfold, modeldim, nimp)) # Initialize array to store imputed validation folds across imputations


      # Compute MI with m=M in the whole dataset, run an imputation with m = M
      imp_data <- mice(combdat,
        m = nimp,
        maxit = 50,
        printFlag = FALSE)

      for (m in 1:nimp) {
        data_compl <- complete(imp_data, m)  # select the complete model

        # Fit a Logistic Cox model on the completed training data, remove the validation data portion
        fit <- glm(formula , family = family, data = data_compl[-(1:lengthfold), ])

        #X_m[m, ] <- model.matrix(formulapredict, data = data_compl[folddef[[k]], ])

        X_m[, , m]  <- model.matrix(formulapredict, data = data_compl[(1:lengthfold),]) # Save the mth completed validation model matrix

        coef_m[m, ] <- fit$coefficients       # Save the coefficients
      }

      coefs[k, ] <-
        apply(coef_m, 2, mean)     # Pool the coefficients (and save them)
      #Xb[folddef[[k]], ] <-
      #  t(X_m %*% coefs[k, ]) # Compute X*b for the validation set data
      Xb[folddef[[k]], ] <- apply((X_m)*((rep(1,lengthfold))%*%t(coefs[k,]))%o%(rep(1,nimp)),c(1,3),sum) # Apply pooled coefficients to all imputed validation data matrices and save linear predictors
    }

    Pred <-
      .expit(Xb) # Compute the probability for each individual in the validation set and save the predictions

    output <- list(pred = Pred, linpred = Xb)
    output
  }
