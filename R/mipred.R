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
#' @param folds Number of fold-partitions defined within \code{newdata}.
#' An integer from 1 to \code{m}. Defaults to NULL which internally
#' sets \code{folds=m}, which puts each observation in \code{newdata}
#' into its own singleton fold. The minimum value \code{folds=1}
#'  would predict the entire set \code{newdata} in a single step without partitioning.
#' @param method Imputation combination method. This defaults to
#'   \code{"averaging"} for the prediction-averaging approach. The alternative
#'   \code{"rubin"} applies the Rubin's rules pooled model.
#' @param mice.options Optional list containing arguments to be supplied to \code{mice}. Refer to the \code{mice} documentation for details.
#' The following options may be specified: \code{method}, \code{predictorMatrix}, \code{blocks},
#' \code{visitSequence}, \code{formulas}, \code{blots}, \code{post}, \code{defaultMethod},
#' \code{maxit}, \code{printFlag}, \code{seed}, \code{data.init}. Please refer to the
#' \code{mice} documentation for the description of these options. To set the number
#' of imputations \code{nimp} should be used. \code{seed} may be specified as a numeric vector
#' of length \code{nimp*folds} when \code{method} is set to \code{averaging} and of length \code{folds}
#' when \code{method} is set to \code{rubin}. Setting \code{seed} to a vector will cause each next
#' call to \code{mice} to use the next seed value in the vector. Setting the seed to a single
#' numeric value will cause all instances of
#' mice to use that same seed value. If you specify a seed vector of insufficient length
#' then the values will be recycled. The required length is \code{folds*nimp} for the averaging
#' approach and length \code{folds} for the rubin approach. The \code{defaultMethod} is set to
#' \code{c("pmm", "logreg", "polyreg", "polr")} by default. The default setting for
#' \code{printFlag} is FALSE. The default for \code{maxit} is 50. All other options are set
#' to \code{NULL} by default.
#'
#' @return A list containing predictions.
#' \describe{ \item{\code{pred}}{Matrix
#'   of predictions on the scale of the response variable of dimension \code{m}
#'   by \code{nimp}.} \item{\code{linpred}}{Matrix of predictions on the scale
#'   of the linear predictor of dimension \code{m} by \code{nimp}.} }
#'
#' @author Bart J A Mertens, \email{b.mertens@lumc.nl}
#' @references \url{https://arxiv.org/abs/1810.05099}
#' @seealso \code{\link{mice}}
#'
#' @examples
#' \dontrun{
#' # Generate a copy of the cll data and construct binary outcome from survival information
#' cll_bin<-cll
#' cll_bin$srv5y_s[cll_bin$srv5y>12] <- 0  # Apply administrative censorship at t=12 months
#' cll_bin$srv5y[cll_bin$srv5y>12]  <- 12
#' cll_bin$Status[cll_bin$srv5y_s==1]<- 1  # Define the new binary "Status" outcome variable
#' cll_bin$Status[cll_bin$srv5y_s==0] <- 0  # As numeric -> 1:Dead, 0:Alive
#' cll_bin$Censor <- NULL # Remove survival outcomes
#' cll_bin$srv5y <- NULL
#' cll_bin$srv5y_s <- NULL
#'
#' # Predict observations 501 to 504 using the first 100 records to calibrate predictors
#' # Remove the identification variable before prediction calibration and imputation.
#' # Remove outcome for new observations
#' # Apply prediction-averaging using 5 imputations, set mice option maxit=5.
#' # Note these settings are only for illustration and should be set to higher values for
#' # practical use, particularly for nimp.
#' output<-mipred(Status ~ age10+cyto, family=binomial, data=cll_bin[1:100,-1],
#'   newdata=cll_bin[501:504,c(-1,-10)], nimp=5, mice.options=list(maxit=5))
#' }
#'
#' @export
#' @importFrom mice mice complete
#' @importFrom stats glm model.matrix
#'
mipred <-
  function(formula,
    family,
    data,
    newdata,
    nimp,
    folds = NULL,
    method = "averaging",
    mice.options = NULL) {
    call <- match.call()

    if (missing(formula))
      stop("formula not provided")
    if (!inherits(formula, "formula"))
      stop("formula argument not recognized as formula")

    if (is.character(family))
      family <-
        get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("family not recognized")
    }

    if (missing(data))
      stop("data not provided")
    if (!is.data.frame(data))
      stop("data should be a data frame", call. = FALSE)
    dup <- duplicated(colnames(data))
    if (any(dup))
      stop("Duplicate names found in data: ",
        paste(colnames(data)[dup], collapse = ", "))

    if (missing(newdata))
      stop("newdata not provided")
    if (!is.data.frame(newdata))
      stop("newdata should be a data frame", call. = FALSE)
    dupnew <- duplicated(colnames(newdata))
    if (any(dupnew))
      stop("Duplicate names found in newdata: ",
        paste(colnames(newdata)[dupnew], collapse = ", "))
    rownames.newdata <- rownames(newdata)

    if (missing(nimp))
      stop("nimp not provided")
    if (!is.numeric(nimp) | !is.vector(nimp)) {
      stop("nimp must be integer number greater than zero")
    } else {
      if (!(length(nimp) == 1)) {
        stop("nimp must be integer number greater than zero")
      } else{
        if ((as.integer(nimp) != nimp) | (nimp < 1))
          stop("nimp must be integer number greater than zero")
      }
    }

    if (is.null(folds)) {
      folds <-
        nrow(newdata) # force default to leave-one-out when NULL
    } else {
      if (!is.numeric(folds) | !is.vector(folds)) {
        stop("folds must be integer number greater than zero")
      } else {
        if (!(length(folds) == 1)) {
          stop("folds must be integer number greater than zero")
        } else{
          if ((as.integer(folds) != folds) | (folds < 1))
            stop("folds must be integer number greater than zero")
        }
      }
    }
    if (folds > nrow(newdata))
      stop("folds too large")

    if (is.null(mice.options)) {
      # setting default mice options
      mice.options[["method"]] <- NULL
      mice.options[["where"]] <- NULL
      mice.options[["visitSequence"]] <- NULL
      mice.options[["blots"]] <- NULL
      mice.options[["post"]] <- NULL
      mice.options[["defaultMethod"]] <-
        c("pmm", "logreg", "polyreg", "polr")
      mice.options[["maxit"]] <- 50
      mice.options[["printFlag"]] <- FALSE
      mice.options[["seed"]] <- NA
      mice.options[["data.init"]] <- NULL
    } else {
      # updating user provided mice options
      if (is.null(mice.options[["defaultMethod"]]))
        mice.options$defaultMethod <-
          c("pmm", "logreg", "polyreg", "polr")
      if (is.null(mice.options[["maxit"]]))
        mice.options$maxit <- 50
      if (is.null(mice.options[["printFlag"]]))
        mice.options$printFlag <- FALSE
      if (is.null(mice.options[["seed"]]))
        mice.options$seed <- NA
    }
    if (!is.null(mice.options[["m"]]))
      warning("mice.options argument m ignored, using nimp instead")

    if (all(!is.na(mice.options[["seed"]]))) {
      if (!is.vector(mice.options[["seed"]])) {
        stop("in mice.options seed must be a vector")
      } else {
        if (!is.numeric(mice.options[["seed"]])) {
          stop("in mice.options seed must be a numeric vector")
        } else {
          if (!all(as.integer(mice.options[["seed"]]) == mice.options[["seed"]]) |
              any(mice.options[["seed"]] < 0)) {
            stop("in mice.options seed must be positive integer")
          } else{
            if (method == "rubin") {
              mice.options[["seed"]] <-
                matrix(mice.options[["seed"]],
                  nrow = folds,
                  ncol = 1)[, 1]
            } else {
              mice.options[["seed"]] <-
                matrix(mice.options[["seed"]],
                  nrow = folds * nimp,
                  ncol = 1)[, 1]
            }
          }
        }
      }
    }

    # generate calibration data.frame
    mf <- match.call(expand.dots = FALSE)
    m <-
      match(c("formula", "data", "subset", "weights"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action = quote(na.pass)
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <-
      quote(stats::model.frame) #mf[[1L]] <- as.name("model.frame") # make a model frame
    mf[[2L]] <-
      as.formula(paste(mf[[2L]][[2L]], "~.")) # keep all predictors
    YX <-
      eval(mf, parent.frame()) # data.frame containing response and ALL predictors

    # generate validation data.frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "newdata"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$na.action = quote(na.pass)
    mf$drop.unused.levels <- TRUE
    names(mf)[3L] <-
      "data" # change to "data" to setup call to model.frame for predictor data
    mf[[2L]] <- as.formula("~.") # keep all predictors
    mf[[1L]] <-
      quote(stats::model.frame) #mf[[1L]] <- as.name("model.frame")
    Xnew <-
      eval(mf, parent.frame()) # data.frame containing ALL predictors only

    YX.vars <- attributes(attributes(YX)$terms)$term.labels
    Xnew.vars <- attributes(attributes(Xnew)$terms)$term.labels

    if (!setequal(YX.vars, Xnew.vars))
      #check for predictors compatibility between calibration and validation data
      stop("data and newdata not compatible: check consistency of included predictors")

    # In this release not yet survival prediction
    if (method == "averaging") {
      # prediction-averaging method
      output <-
        .glm.mipred.cmb1(formula, family, YX, Xnew, nimp, folds, mice.options)
    } else {
      # Rubin's rules model
      output <-
        .glm.mipred.cmb2(formula, family, YX, Xnew, nimp, folds, mice.options)
    }

    row.names(output$pred) <-
      row.names(output$linpred) <- rownames.newdata
    output <- list(call = call,
      pred = output$pred,
      linpred = output$linpred)
    output
  }
