#' Expit function converting odds to probability
#'
#' @param x Probability vector

.expit <- function(x){exp(x)/(1+exp(x))}

#' Binary prediction using multiple imputation - prediction-averaging method
#'
#' @param formula Formula used by fitting and prediction method
#' @param family Error distribution also determining the link function used
#' @param dataset A data frame containing calibration data
#' @param newdata A dataframe containing observations to be predicted
#' @param nimp Number of imputations for each observation

.glm_mipred_cmb1 <- function(formula, family, dataset, newdata, nimp){

  formulapredict<-formula
  formulapredict[[2]]<-NULL # use delete.response?
  newn <- nrow(newdata)

  newdata<-base::cbind(NA,newdata) # append column vector of NA for absent response in validation set
  names(newdata)[[1]]<-names(dataset)[[1]] # rename to allow rbind below

  # Define some matrices we need later
  Xb  <- matrix(NA, nrow = newn, ncol = nimp)   # linear predictor for new observations
  Pred <- matrix(NA, nrow = newn, ncol = nimp)  # predicted probabilities for new observations

  for (k in 1:newn) {

    combdat <- base::rbind(newdata[k, ], dataset)  # case-by-case predictions for new observations

    for (m in 1:nimp) {

      # Compute MI (once) in the combined dataset, run an imputation with m = 1
      imp_data <- mice(combdat, m = 1, maxit = 50, print = FALSE)
      data_compl <- complete(imp_data, 1)  # select the completed model data

      # Fit a Logistic Cox model on the completed training data, remove the validation data portion
      fit  <- glm(formula, family = family, data = data_compl[-1, ]) # default family=binomial("logit")
      coef <- fit$coefficients # Save the model coefficients

      X_m <- model.matrix(formulapredict, data = data_compl[1,])   # Define the regression matrix and remove status next
      Xb[k,m] <- X_m %*% coef            # Compute X*b only for the validation sample
    }
#    print(k)
  }

  # Compute the probability for each individual in the validation set
  Pred <- .expit(Xb)

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

.glm_mipred_cmb2 <- function(formula, family, dataset, newdata, nimp){

  formulapredict<-formula
  formulapredict[[2]]<-NULL   # use delete.response?
  newn <- nrow(newdata)

  modeldim<-dim(model.matrix(formula , data=dataset))[2] # length of model coefficient vector
  newdata<-base::cbind(NA,newdata) # append column vector of NA for absent response in validation set
  names(newdata)[[1]]<-names(dataset)[[1]] # rename to allow rbind below

  # Define some matrices we need later
  coef_m <- matrix(NA, nrow=nimp, ncol=modeldim)
  X_m <- matrix(NA, nrow=nimp, ncol=modeldim)
  Xb <- matrix(NA, nrow=newn, ncol=nimp)
  Pred <- matrix(NA, nrow=newn, ncol=nimp)

  coefs  <- matrix(NA, nrow=newn, ncol=modeldim) # only if we want to save the results

  for (k in 1:newn) {

    combdat <- base::rbind(newdata[k, ], dataset)  # case-by-case predictions for new observations

    # Compute MI with m=M in the whole dataset, run an imputation with m = M
    imp_data <- mice(combdat, m = nimp, maxit=50, print=FALSE)

    for (m in 1:nimp){

      data_compl <- complete(imp_data, m)  # select the complete model

      # Fit a Logistic Cox model on the completed training data, remove the validation data portion
      fit <- glm(formula , family=family, data=data_compl[-1,])

      X_m[m,] <- model.matrix(formulapredict, data=data_compl[1,])

      coef_m[m,] <- fit$coefficients       # Save the coefficients
    }

    coefs[k,] <- apply(coef_m, 2, mean)     # Pool the coefficients (and save them)
    Xb[k,] <-t(X_m%*%coefs[k,]) # Compute X*b for the validation set data
#    print(k)
  }

  Pred <- .expit(Xb) # Compute the probability for each individual in the validation set and save the predictions

  output <- list(pred = Pred, linpred = Xb)
  output
}






