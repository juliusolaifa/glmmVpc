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

