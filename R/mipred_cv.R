#' Cross-validation prediction using multiple imputation
#'
#' Calculates cross-validated predictions based on within-sample  assessment and calibration using
#' generalized linear models with multiple imputations to account for missing values
#' in predictor data.
#'
#' @param formula A formula object providing a symbolic description of the
#'   prediction model to be fitted.
#' @param family Specification of an appropriate error distribution and link
#'   function.
#' @param data A data.frame containing calibration data on \code{n} samples.
#'   Variables declared in \code{formula} must be found in \code{data}.
#' @param nimp Number of imputations used in the prediction of each observation.
#' @param folds Number of fold-partitions used in cross-validation. Put \code{folds=n}
#' for leave-one-out cross-validation.
#' @param method Imputation combination method. This defaults to
#'   \code{"averaging"} for the prediction-averaging approach. The alternative
#'   \code{"rubin"} applies the Rubin's rules pooled model.
#'
#' @return A list containing predictions.
#' \describe{\item{\code{pred}}{Matrix of predictions on the scale of the response
#' variable of dimension \code{n} by \code{nimp}.}
#' \item{\code{linpred}}{Matrix of predictions on the scale
#'   of the linear predictor of dimension \code{n} by \code{nimp}.}}
#'
#' @author Bart J A Mertens, \email{b.mertens@lumc.nl}
#' @references \url{https://arxiv.org/abs/1810.05099}
#' @seealso \code{\link{mice}}
#'
#' @examples
#' cll_bin<-cll # Generate a copy and construct binary outcome from survival information
#' cll_bin$srv5y_s[cll_bin$srv5y>12] <- 0  # Apply administrative censorship at t=12 months
#' cll_bin$srv5y[cll_bin$srv5y>12]  <- 12
#' cll_bin$Status[cll_bin$srv5y_s==1]<- 1  # Define the new binary "Status" outcome variable
#' cll_bin$Status[cll_bin$srv5y_s==0] <- 0  # As numeric -> 1:Dead, 0:Alive
#' cll_bin$Censor <- NULL # Remove survival outcomes
#' cll_bin$srv5y <- NULL
#' cll_bin$srv5y_s <- NULL
#'
#' # Cross-validate prediction using logistic regression in the first 100 samples
#' # Apply prediction-averaging using 5 imputations and 5 folds
#' output<-mipred.cv(Status~perfstat+remstat+cyto, family=binomial, data=cll_bin[1:100,],
#' nimp=5, folds=5)
#'
#' @export
#' @importFrom mice mice complete
#' @importFrom stats glm model.matrix
#'
mipred.cv <- function(formula, family, data, nimp, folds, method = "averaging"){

  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)

  # generate calibration data.frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action=quote(na.pass)
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame) #mf[[1L]] <- as.name("model.frame")
  YX <- eval(mf, parent.frame()) # data.frame containing response and predictors

  # In this preliminary release only binary prediction
  if (method=="averaging"){ # prediction-averaging method
    output <- .glm_mipred_cmb1_cv(formula, family, YX, nimp, folds)
  } else { # Rubin's rules model
    output <- .glm_mipred_cmb2_cv(formula, family, YX, nimp, folds)
  }
  output<-list(call=call,pred=output$pred,linpred=output$linpred)
  output
}


