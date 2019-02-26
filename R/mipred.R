#' Prediction using multiple imputation
#'
#' Calculates predictions from generalized linear models when multiple
#' imputations are used to account for missing values in predictor data.
#'
#' @param formula A formula object providing a symbolic description of the
#'   prediction model to be fitted.
#' @param family Specification of an appropriate error distribution and link
#'   function.
#' @param data A data.frame containing calibration data on \code{n} samples.
#'   Variables declared in \code{formula} must be found in \code{data}.
#' @param newdata A data.frame containing the predictors for observations to be
#'   predicted on \code{m} samples. This must have the same structure and
#'   variables as \code{data}, except for the outcome variable which is ignored
#'   in the construction of the predictions and can therefor be excluded from
#'   the object.
#' @param nimp Number of imputations used in the prediction of each observation.
#' @param method Imputation combination method. This defaults to
#'   \code{"averaging"} for the prediction-averaging approach. The alternative
#'   \code{"rubin"} applies the Rubin's rules pooled model.
#'
#' @return A list containing predictions. \describe{ \item{\code{pred}}{Matrix
#'   of predictions on the scale of the response variable of dimension \code{m}
#'   by \code{nimp}.} \item{\code{linpred}}{Matrix of predictions on the scale
#'   of the linear predictor of dimension \code{m} by \code{nimp}.} }
#'
#' @author Bart J A Mertens, \email{b.mertens@lumc.nl}
#' @references \url{https://arxiv.org/abs/1810.05099}
#' @seealso \code{\link{mice}}
#'
#' @examples
#'
#' #load("cll.rda")
#' cll_bin<-cll # Generate a copy and construct binary outcome from survival information
#' cll_bin$srv5y_s[cll_bin$srv5y>12] <- 0  # Apply administrative censorship at t=12 months
#' cll_bin$srv5y[cll_bin$srv5y>12]  <- 12
#' cll_bin$Status[cll_bin$srv5y_s==1]<- 1  # Define the new binary "Status" outcome variable
#' cll_bin$Status[cll_bin$srv5y_s==0] <- 0  # As numeric -> 1:Dead, 0:Alive
#' cll_bin$Censor <- NULL # Remove survival outcomes
#' cll_bin$srv5y <- NULL
#' cll_bin$srv5y_s <- NULL
#'
#' set.seed(12345)
#' # Predict observations 501 to 507 using the first 100 records to calibrate predictors
#' # Apply prediction-averaging using 5 imputations
#' output<-mipred(Status~perfstat+remstat+cyto, family=binomial, data=cll_bin[1:100,],
#'   newdata=cll_bin[501:507,-10], nimp=5)
#'
#' @export
#' @importFrom mice mice complete
#' @importFrom stats glm model.matrix
#'
mipred <- function(formula, family, data, newdata, nimp, method="averaging")  {

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

  # generate validation data.frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "newdata"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$na.action=quote(na.pass)
  mf$drop.unused.levels <- TRUE
  names(mf)[3]<-"data" # change to "data" to setup call to model.frame for predictor data
  mf[[2]][[2]]<-NULL # remove response from formula to generate data.frame with predictors only
  mf[[1L]] <- quote(stats::model.frame) #mf[[1L]] <- as.name("model.frame")
  Xnew <- eval(mf, parent.frame()) # contains predictors only

  # In this preliminary release only binary prediction
  if (method=="averaging"){ # prediction-averaging method
    output <- .glm_mipred_cmb1(formula, family, YX, Xnew, nimp)
  } else { # Rubin's rules model
    output <- .glm_mipred_cmb2(formula, family, YX, Xnew, nimp)
  }
  output<-list(call=call,pred=output$pred,linpred=output$linpred)
  output
}


