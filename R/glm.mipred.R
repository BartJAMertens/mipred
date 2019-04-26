#' Generalized linear model prediction using multiple imputation - prediction-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param newdata A dataframe containing observations to be predicted
#' @param nimp Number of imputations for each observation
#' @param folds Number of folds defined in newdata
#' @param miop Mice options

.glm.mipred.cmb1 <-
  function(formula,
    family,
    dataset,
    newdata,
    nimp,
    folds,
    miop) {
    newn <- nrow(newdata)

    newdata <-
      base::cbind(NA, newdata) # append column vector of NA for absent response in validation set
    names(newdata)[[1]] <-
      names(dataset)[[1]] # rename to allow rbind below

    # Define some matrices we need later
    Xb <-
      Pred <-
      matrix(NA, nrow = newn, ncol = nimp)   # Initialize predictor matrices for new observations

    for (m in 1:nimp) {
      suppressWarnings(folddef <- split(sample(newn, newn, replace = F), as.factor(1:folds)))

      for (k in 1:folds) {
        folddef[[k]] <-
          sort(folddef[[k]])    # Sort the ids of the k-th fold
        combdat <-
          base::rbind(newdata[folddef[[k]],], dataset)  # case-by-case predictions for new observations
        lengthfold <- length(folddef[[k]])

        imp_data <-
          impute(combdat, miop, 1, miop[["seed"]][k + (m - 1) * folds])
        data_compl <-
          complete(imp_data, 1)  # select the completed model data

        # Fit a model on the completed training data, remove the validation data portion
        fit  <-
          glm(formula, family = family, data = data_compl[-(1:lengthfold),])
        coefs <- fit$coefficients # Save the model coefficients

        X_m <-
          model.matrix(formula,
            data = data_compl[(1:lengthfold), ],
            drop.unused.levels = FALSE) # Define the regression matrix and remove status next

        Xb[folddef[[k]], m] <-
          X_m[, names(fit$coefficients)] %*% coefs  # Compute X*b only for the validation sample

        Pred[folddef[[k]], m] <-
          predict.glm(fit, data_compl[(1:lengthfold), ], type = "response") # This will generate error if more factor levels than calibration sets
      }
      # print(m)  change this to verbose option
    }

    output <- list(pred = Pred, linpred = Xb)
    output
  }


#' Generalized linear model prediction using multiple imputation - Rubin's rule coefficient-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param newdata A dataframe containing observations to be predicted
#' @param nimp Number of imputations for each observation
#' @param folds Number of folds defined in newdata
#' @param miop Mice options

.glm.mipred.cmb2 <-
  function(formula,
    family,
    dataset,
    newdata,
    nimp,
    folds,
    miop) {
    attributes(dataset)$terms<-NULL
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
    Xb <-
      Pred <-
      matrix(NA, nrow = newn, ncol = nimp) # Initialize predictor matrices for new observations

    coefs  <-
      matrix(NA, nrow = folds, ncol = modeldim) # only if we want to save the results
    suppressWarnings(folddef <- split(sample(newn, newn, replace = F), as.factor(1:folds)))

    for (k in 1:folds) {
      folddef[[k]] <-
        sort(folddef[[k]])    # Sort the ids of the k-th fold
      lengthfold <- length(folddef[[k]])
      combdat <-
        base::rbind(newdata[folddef[[k]],], dataset)  # case-by-case predictions for new observations
      X_m <-
        array(NA, dim = c(lengthfold, modeldim, nimp)) # Initialize array to store imputed validation folds across imputations

      # Compute MI with m=nimp in the whole dataset, run an imputation with m = nimp, use previously generated options
      imp_data <- impute(combdat, miop, nimp, miop[["seed"]][k])
# print(imp_data$loggedEvents) change to Verbose option

      for (m in 1:nimp) {
        data_compl <- complete(imp_data, m)  # select the complete model

        # Fit a Logistic Cox model on the completed training data, remove the validation data portion
        fit <-
          glm(formula , family = family, data = data_compl[-(1:lengthfold), ])

        X_m[, , m]  <-
          model.matrix(formula,
            data = data_compl[(1:lengthfold),],
            drop.unused.levels = FALSE) # Save the mth completed validation model matrix
        coef_m[m, ] <-
          fit$coefficients # Save the coefficients
      }

      coefs[k, ] <-
        apply(coef_m, 2, mean)     # Pool the coefficients (and save them)
      Xb[folddef[[k]], ] <-
        apply((X_m) * ((rep(1, lengthfold)) %*% t(coefs[k,])) %o% (rep(1, nimp)), c(1, 3), sum) # Apply pooled coefficients to all imputed validation data matrices and save linear predictors

      fit$coefficients <-
        coefs[k, ] # replace coefficient vector with pooled coefficients

      if (!any(is.na(newdata[folddef[[k]], -1]))) {
        Pred[folddef[[k]], ] <-
          matrix(predict.glm(fit, data_compl[(1:lengthfold),], type = "response"), ncol =
              1) %*% rep(1, nimp)
      } else {
        for (m in 1:nimp) {
          data_compl <- complete(imp_data, m)
          Pred[folddef[[k]], m] <-
            predict.glm(fit, data_compl[(1:lengthfold),], type = "response") # This will generate error if more factor levels than for calibration sets
        }
      }
      # print(k) change this to Verbose option
    }
    output <- list(pred = Pred, linpred = Xb)
    output
  }


#' General imputation routine for mipred
#'
#' @param combdat Dataset to be imputed
#' @param miop Mice options list
#' @param nimp Number of imputations
#' @param seed Single numerical seed value

impute <- function(combdat, miop, nimp, seed) {
  mioppredictorMatrix <- is.null(miop$predictorMatrix)
  miopblocks <- is.null(miop$blocks)
  miopformulas <- is.null(miop$formulas)

  # check predictorMatrix, blocks and formulas
  if (mioppredictorMatrix & miopblocks & miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        #predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        #blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        #formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (!mioppredictorMatrix &
      miopblocks & miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        #blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        #formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (mioppredictorMatrix & !miopblocks & miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        #predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        #formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (mioppredictorMatrix & miopblocks & !miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        #predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        #blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (!mioppredictorMatrix &
      !miopblocks & miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        #formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (!mioppredictorMatrix &
      miopblocks & !miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        #blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (mioppredictorMatrix &
      !miopblocks & !miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        #predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  if (!mioppredictorMatrix &
      !miopblocks & !miopformulas) {
    imp_data <-
      mice(
        combdat,
        m = nimp,
        method = miop[["method"]],
        predictorMatrix = miop[["predictorMatrix"]],
        where = miop[["where"]],
        blocks = miop[["blocks"]],
        visitSequence =  miop[["visitSequence"]],
        formulas = miop[["formulas"]],
        blots = miop[["blots"]],
        post = miop[["post"]],
        defaultMethod = miop[["defaultMethod"]],
        maxit = miop[["maxit"]],
        printFlag = miop[["printFlag"]],
        seed = seed,
        data.init = miop[["data.init"]]
      )
  }
  imp_data
}


#' Expit function converting odds to probability
#'
#' @param x Probability vector

.expit <- function(x) {
  exp(x) / (1 + exp(x))
}
