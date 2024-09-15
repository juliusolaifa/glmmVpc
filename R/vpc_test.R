vpc_test <- function(fitObj, null_formula = NULL,  test_type = c("LRT", "Score")) {
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
