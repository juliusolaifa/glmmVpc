#' @export
#' @method vpc.test glmmfit
#' @rdname vpc.test
vpc.test.glmmfit <- function(fitObj, null_formula, type=c("classical",
                                                          "self",
                                                          "zhang",
                                                          "julius",
                                                          "all"), ...) {

  type <- match.arg(type)
  if(inherits(fitObj, "vpcObj")) {
    fitObj <- fitObj$modObj
  }
  data <- stats::model.frame(fitObj)
  family <- fitObj$family

  fit0 <- singleGLMMFit(formula = null_formula, data = data, family = family)

  logLik0 <- stats::logLik(fit0) #-fit0$modObj$fit$objective#
  logLik1 <- stats::logLik(fitObj) #-fitObj$modObj$fit$objective#

  test_stat <- 2 * (as.numeric(logLik1) - as.numeric(logLik0))

  if(type == "classical") {
    df_diff <- attr(logLik1, "df") - attr(logLik0, "df")
    if (df_diff <= 0) {
      stop("Null model is not nested within the alternative model.")
    }
    p_value <- stats::pchisq(test_stat, df = df_diff, lower.tail = FALSE)
  } else {
    config <- obtain_config(fit0, fitObj)
    p <- pr11(fitObj, type=type)
    if(type == "all") {
      df_diff <- attr(logLik1, "df") - attr(logLik0, "df")
      types <- c("self", "zhang", "julius")
      p_values <- lapply(types, function(t) adj_chisq(test_stat = test_stat,
                                                      config = config, type = t,
                                                      p = p))
      return(list(LR_stat = test_stat,
                  p_value.c = stats::pchisq(test_stat, df = df_diff, lower.tail = FALSE),
                  p_value.s = p_values[[1]], p_value.z = p_values[[2]], p_value.j = p_values[[3]]))
    }else{
      p_value <- adj_chisq(test_stat = test_stat, config=config, type=type, p=p)
    }
  }

    return(list(LR_stat = test_stat, p_value = p_value))

}

#' @export
#' @method vpc.test vpcObj
vpc.test.vpcObj <- function(fitObj, null_formula, type=c("classical",
                                                         "self",
                                                         "zhang",
                                                         "julius",
                                                         "all"), ...) {
  vpc.test.glmmfit(fitObj, null_formula, type, ...)
}


#' @export
#' @method vpc.test Glmmfits
#' @rdname vpc.test
vpc.test.Glmmfits <- function(fitObj, null_formula,
                              type=c("classical",
                                     "self",
                                     "zhang",
                                     "julius",
                                     "all"), num_cores = 1, ...) {
  n <- length(fitObj)
  type <- match.arg(type)
  results <- data.frame(t(pbapply::pbsapply(seq_along(fitObj), function(i) {

    lhs <- paste0(deparse(null_formula[[2]]), i)
    rhs <- deparse(null_formula[[3]])
    formula_string <- paste(lhs, "~", rhs)
    dynamic_formula <- stats::as.formula(formula_string)
    # formula_string <- paste(paste0(deparse(null_formula[[2]]),i), "~", deparse(null_formula[[3]]))
    # dynamic_formula <- stats::as.formula(formula_string)
    vpc.test(fitObj = fitObj[[i]], null_formula = dynamic_formula, type=type)
  }, cl = num_cores)))
  #results$p_value <- as.numeric(results$p_value)
  #results$LR_stat <- as.numeric(results$LR_stat)
  return(results)
}

#' @export
#' @method vpc.test VpcObj
vpc.test.VpcObj <- function(fitObj, null_formula, type=c("classical",
                                                         "self",
                                                         "zhang",
                                                         "julius",
                                                         "all"), ...) {
  vpc.test.Glmmfits(fitObj, null_formula, type, ...)
}


#' Perform Variance Partition Coefficient (VPC) Test
#'
#' This function calculates the Variance Partition Coefficient (VPC) for a fitted model object.
#' The specific method used will depend on the class of the provided model.
#'
#' @param fitObj The fitted model object of class \code{glmmfit} or a list of fitted model objects
#' of class \code{glmmfits} for which the VPC is to be calculated.
#' @param null_formula A formula to specify a null model for comparison.
#' @param type Character string indicating the method used to calculate the chi-square distribution.
#' Options are \code{"classical"}, \code{"self"}, \code{"zhang"}, \code{"julius"}, or \code{"all"}.
#' @param num_cores (numeric) Number of cores to use for parallel processing.
#' Only applicable when \code{fitObj} is a list of fitted models (i.e., class \code{glmmfits}).
#' @param ... Additional arguments passed to specific methods.
#'
#' @return A numeric value or a list representing the calculated VPC for the given model,
#' including additional details if applicable.
#' @export
vpc.test <- function(fitObj, null_formula, type=c("classical",
                                                  "self",
                                                  "zhang",
                                                  "julius",
                                                  "all"), ...) {
  UseMethod("vpc.test")
}


#' Adjusted Chi-Square Test Statistic
#'
#' This function calculates the p-value for a chi-square test statistic based on the
#' configuration parameters passed in `config`. It supports different cases for
#' parameters on and off the boundary.
#'
#' @param test_stat Numeric value representing the chi-square test statistic.
#' @param config A list containing the configuration parameters:
#' @param ... additional parameters, like `p` probability of being in the admissible
#' region. See \url{https://www.jstor.org/stable/2289471}
#'
#' \describe{
#'   \item{q}{Parameter of interest on the boundary.}
#'   \item{r}{Parameter of interest not on the boundary.}
#'   \item{s}{Nuisance parameter on the boundary (in this function, it is always 0).}
#'   \item{t}{Nuisance parameter not on the boundary.}
#' }
#'
#' @return A numeric value representing the calculated p-value based on the
#' specified configuration.
#'
#' @details
#' The function evaluates the test statistic based on three cases:
#' \describe{
#'   \item{Case I}{If \code{q == 0} and \code{s == 0}, it calculates the p-value
#'   using \code{r} degrees of freedom.}
#'   \item{Case II}{If \code{q == 1}, \code{r == 0}, and \code{s == 0}, it
#'   calculates a weighted p-value combining \code{df = 1} and \code{test_stat == 0}.}
#'   \item{Case III}{If \code{q == 1} and \code{s == 0}, it calculates a weighted
#'   p-value using degrees of freedom \code{r - 1} and \code{r}.}
#' }
#'
#' @references
#' Self, S.G., & Liang, K. (1987). *On Asymptotic properties of maximum
#' likelihood estimators
#' and likelihood ratio tests under non-standard conditions*.
#' Journal of American Statistical Association, 82(398), 605-610-145.
#' \url{https://www.jstor.org/stable/2289471}
#'
#' @export
adj_chisq <- function(test_stat, config, ...) {

  # Ensure config has all the necessary components
  if (is.null(config$q) || is.null(config$r) || is.null(config$s) || is.null(config$t)) {
    stop("Config must contain 'q', 'r', 's', and 't'.")
  }

  calculate_p_value <- function(q, r, s, test_stat) {
    caseI <- q == 0 && s == 0
    caseII <- q == 1 && r == 0 && s == 0
    caseIII <- q == 1 && r == 1 && s == 0
    caseIV <- q == 2


    if (caseI)   {
      return(stats::pchisq(test_stat, df = r, lower.tail = FALSE))
    } else if (caseII) {
      chi0 <- test_stat == 0
      chi1 <- stats::pchisq(test_stat, df = 1, lower.tail = FALSE)
      return(0.5 * chi0 + 0.5 * chi1)
    } else if (caseIII) {
      chi1 <- stats::pchisq(test_stat, df = 1, lower.tail = FALSE)
      chi2 <- stats::pchisq(test_stat, df = 2, lower.tail = FALSE)
      return(0.5 * chi1 + 0.5 * chi2)
    } else if(caseIV) {
      args <- list(...)
      p <- args$p
      type <- args$type
      chi0 <- test_stat == 0
      chi1 <- stats::pchisq(test_stat, df = 1, lower.tail = FALSE)
      chi2 <- stats::pchisq(test_stat, df = 2, lower.tail = FALSE)
      chi3 <- stats::pchisq(test_stat, df = 3, lower.tail = FALSE)
      pval.sz <- (0.5-p) * chi0 + 0.5 * chi1 + p*chi2
      pval.j <- (0.5-p) * chi1 + 0.5 * chi2 + p*chi3
      if (type == "self" || type == "zhang") {
        return(pval.sz)
      } else if(type == "julius"){
        return(pval.j)
      }
    }else {
      stop("Invalid configuration for q, r, s.")
    }
  }

  p_value <- calculate_p_value(config$q, config$r, config$s, test_stat)

  return(p_value)
}



#' Obtain Configuration Parameters for Boundary and Non-Boundary Models
#'
#' This function calculates the parameters `q`, `r`, `s`, and `t` based on the
#' difference
#' between a null and alternative model, particularly in cases where certain
#' parameters
#' are on the boundary or not.
#'
#' @param null A model object representing the null hypothesis. This object
#' should contain
#'   `modObj` and `Sigma`, where `modObj` is the model object and `Sigma` is a
#'   matrix indicating
#'   the parameters on the boundary.
#' @param alt A model object representing the alternative hypothesis. This
#' object should contain
#'   `modObj` and `Sigma`, where `modObj` is the model object and `Sigma` is a
#'   matrix indicating
#'   the parameters on the boundary.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{q}{Number of parameters of interest that are on the boundary in the
#'   alternative model but not in the null model.}
#'   \item{r}{Number of parameters of interest that are not on the boundary in
#'   either the null or alternative model.}
#'   \item{s}{Number of nuisance parameters on the boundary (always 0 in these
#'   cases).}
#'   \item{t}{Number of nuisance parameters not on the boundary.}
#' }
#'
#' @references
#' Self, S.G., & Liang, K. (1987). *On Asymptotic properties of maximum
#' likelihood estimators
#' and likelihood ratio tests under non-standard conditions*.
#' Journal of American Statistical Association, 82(398), 605-610-145.
#' \url{https://www.jstor.org/stable/2289471}
#'
#' @export
obtain_config <- function(null, alt){

  # q : par of interest on the boundary
  # r : par of interest not on boundary
  # s : nuisance par on the boundary (for our use cases this is always `0`)
  # t : nuisance par not on the boundary
  p <- full_par.alt <- nrow(stats::vcov(alt$modObj, full=TRUE))
  full_par.null <- nrow(stats::vcov(null$modObj, full=TRUE))

  bounded_par.alt <- nrow(alt$Sigma)
  bounded_par.null <- ifelse(inherits(null$Sigma, "matrix"), nrow(null$Sigma),
                             length(null$Sigma))

  non_bounded.alt  <- full_par.alt - bounded_par.alt
  non_bounded.null <- full_par.null - bounded_par.null

  q <- bounded_par.alt - bounded_par.null
  r <- non_bounded.alt - non_bounded.null
  s <- 0
  t <- p-q-r-s
  return(list("q"=q,"r"=r,"s"=s,"t"=t))
}



