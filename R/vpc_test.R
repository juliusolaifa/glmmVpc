#' @export
#' @method vpc.test glmmfit
#' @rdname vpc.test
vpc.test.glmmfit <- function(fitObj, null_formula = NULL,  test_type = c("LRT", "Score")) {
  test_type <- match.arg(test_type)

  data <- stats::model.frame(fitObj)
  family <- fitObj$family

  if (is.null(null_formula)) {
    null_formula <- stats::as.formula("response ~ 1")
  }

  null_fit <- singleGLMMFit(formula = null_formula, data = data, family = family)

  if (test_type == "LRT") {
    logLik_null <- stats::logLik(null_fit)
    logLik_alt <- stats::logLik(fitObj)

    df_diff <- attr(logLik_alt, "df") - attr(logLik_null, "df")

    if (df_diff <= 0) {
      stop("Null model is not nested within the alternative model.")
    }

    test_stat <- 2 * (as.numeric(logLik_alt) - as.numeric(logLik_null))
    p_value <- stats::pchisq(test_stat, df = df_diff, lower.tail = FALSE)

    return(list(method = "LRT", LR_stat = test_stat, df = df_diff, p_value = p_value))

    } else if (test_type == "Score") {
      score_vec <- 0#glmmTMB::gradients(null_fit$modObj) ## not implemented
      var_cov_matrix <- stats::vcov(null_fit$modObj)

      # (since Σ = I^(-1))
      test_stat <- t(score_vec) %*% var_cov_matrix %*% score_vec
      df <- length(score_vec)

      p_value <- stats::pchisq(test_stat, df = df, lower.tail = FALSE)

      return(list(method = "Score Test", Score_stat = test_stat, df = df, p_value = p_value))
    }
}


#' @export
#' @method vpc.test Glmmfits
#' @rdname vpc.test
vpc.test.Glmmfits <- function(fitObj, null_formula = NULL, test_type = c("LRT", "Score")) {
  n <- length(fitObj)
  test_type <- match.arg(test_type)

  results <- data.frame(t(pbapply::pbsapply(seq_along(fitObj), function(i) {
    formula_string <- paste(paste0(deparse(null_formula[[2]]),i), "~", deparse(null_formula[[3]]))
    dynamic_formula <- stats::as.formula(formula_string)
    vpc.test(fitObj = fitObj[[i]], null_formula = dynamic_formula, test_type = test_type)
  })))

  return(results)
}


#' Perform Variance Partition Coefficient (VPC) Test
#'
#' This function calculates the Variance Partition Coefficient (VPC) for a fitted model object.
#' The specific method used will depend on the class of the provided model.
#'
#' @param fitObj The fitted model object of class `glmmfit` or a list of fitted model objects
#' of class `glmmfits` for which the VPC is to be calculated.
#' @param null_formula An optional formula to specify a null model for comparison.
#' This can be useful for certain statistical tests (e.g., likelihood ratio tests).
#' @param test_type A character vector indicating the type of test to perform. Options are
#' "LRT" for Likelihood Ratio Test or "Score" for Score Test. The default is "LRT".
#'
#' @return A numeric value or a list representing the calculated VPC for the given model,
#' including additional details if applicable.
#' @export
vpc.test <- function(fitObj, null_formula = NULL, test_type = c("LRT", "Score")) {
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
adj_chisq <- function(test_stat, config) {

  # Ensure config has all the necessary components
  if (is.null(config$q) || is.null(config$r) || is.null(config$s) || is.null(config$t)) {
    stop("Config must contain 'q', 'r', 's', and 't'.")
  }

  calculate_p_value <- function(q, r, s, test_stat) {
    caseI <- q == 0 && s == 0
    caseII <- q == 1 && r == 0 && s == 0
    caseIII <- q == 1 && s == 0

    if (caseI)   {
      return(stats::pchisq(test_stat, df = r, lower.tail = FALSE))
    } else if (caseII) {
      return(0.5 * stats::pchisq(test_stat, df = 1, lower.tail = FALSE) +
               0.5 * (test_stat == 0))
    } else if (caseIII) {
      df1 <- r - 1
      df2 <- r
      return(0.5 * stats::pchisq(test_stat, df = df1, lower.tail = FALSE) +
               0.5 * stats::pchisq(test_stat, df = df2, lower.tail = FALSE))
    } else {
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
  bounded_par.null <- nrow(null$Sigma)

  non_bounded.alt  <- full_par.alt - bounded_par.alt
  non_bounded.null <- full_par.null - bounded_par.null

  q <- bounded_par.alt - bounded_par.null
  r <- non_bounded.alt - non_bounded.null
  s <- 0
  t <- p-q-r-s
  return(list("q"=q,"r"=r,"s"=s,"t"=t))
}



