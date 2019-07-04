#' Cross-validation prediction using multiple imputation
#'
#' Calculates cross-validated predictions based on within-sample assessment and calibration using
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
#' @param folds Number of fold-partitions defined within \code{data} used in cross-validation.
#' An integer from 2 to \code{n}. Defaults to NULL which internally
#' sets \code{folds=n}, which puts each observation in \code{data}
#' into its own singleton fold for leave-one-out cross-validation.
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
#' \donttest{
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
#' # Cross-validate prediction using logistic regression in the first 100 samples
#' # Apply prediction-averaging using 5 imputations, 5 folds and maxit=5.
#' # Note these settings are only for illustration and should be set to higher values for
#' # practical use, particularly for nimp.
#' output<-mipred.cv(Status ~ age10+cyto, family=binomial, data=cll_bin[1:100,-1],
#' nimp=5, folds=5, mice.options=list(maxit=5))
#' }
#'
#' @export
#' @importFrom mice mice complete
#' @importFrom stats glm model.matrix as.formula predict.glm
#'
mipred.cv <-
  function(formula,
    family,
    data,
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
    rownames.data <- rownames(data)
    if (missing(nimp))
      stop("nimp not provided")
    if (!is.numeric(nimp) | !is.vector(nimp)) {
      stop("nimp must be integer number greater than zero")
    } else {
      if (!(length(nimp) == 1)) {
        stop("nimp must be integer number greater than zero")
      } else {
        if ((as.integer(nimp) != nimp) | (nimp < 1))
          stop("nimp must be integer number greater than zero")
      }
    }

    if (is.null(folds)) {
      folds <-
        nrow(data) # force default to leave-one-out when NULL
    } else {
      if (!is.numeric(folds) | !is.vector(folds)) {
        stop("folds must be integer number greater than 1")
      } else {
        if (!(length(folds) == 1)) {
          stop("folds must be integer number greater than 1")
        } else{
          if ((as.integer(folds) != folds) | (folds < 2))
            stop("folds must be integer number greater than 1")
        }
      }
    }
    if (folds > nrow(data))
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
          } else {
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

    # In this release not yet survival prediction
    if (method == "averaging") {
      # prediction-averaging method
      output <-
        .glm.mipred.cmb1.cv(formula, family, YX, nimp, folds, mice.options)
    } else {
      # Rubin's rules model
      output <-
        .glm.mipred.cmb2.cv(formula, family, YX, nimp, folds, mice.options)
    }

    row.names(output$pred) <-
      row.names(output$linpred) <- rownames.data
    output <- list(call = call,
      pred = output$pred,
      linpred = output$linpred)
    output
  }
