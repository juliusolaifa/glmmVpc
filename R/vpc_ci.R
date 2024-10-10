gradient_vpc_engine <- function(beta, Sigma, phi, family, x=NULL, power = NULL, method) {

  if (family == "tweedie" && is.null(power)) {
    stop("Parameter 'power' is required for Tweedie family.")
  }

  # Construct the function string dynamically
  beta_str <- construct_beta_param_str(beta)
  Sigma_str <- construct_sigma_param_str(Sigma)
  param_str <- switch(family,
                      "negative_binomial" = paste(beta_str, "phi", Sigma_str, sep = ", "),
                      "tweedie" = paste(beta_str, "phi", Sigma_str, "power", sep = ", ")
  )

  param_values <- switch(family,
                         "negative_binomial" = c(
                           stats::setNames(beta, paste0("b", seq_along(beta) - 1)),
                           phi = phi,
                           mat2vec(Sigma)),
                         "tweedie" = c(
                           stats::setNames(beta, paste0("b", seq_along(beta) - 1)),
                           phi = phi,
                           mat2vec(Sigma),
                           power = power)
  )

  function_string <- switch(family,
                            "negative_binomial" = paste0(
                              "function(", param_str, ") {",
                              "beta_vals <- c(", beta_str, "); ",
                              "sigma_vals <- c(", Sigma_str, "); ",
                              "mean <- mu(beta_vals, x); ",
                              "variance <- sig(sigma_vals,x); ",
                              "vpc.nb(mean, variance, phi); ",
                              "}"
                            ),
                            "tweedie" = paste0(
                              "function(", param_str, ") {",
                              "beta_vals <- c(", beta_str, "); ",
                              "sigma_vals <- c(", Sigma_str, "); ",
                              "mean <- mu(beta_vals, x); ",
                              "variance <- sig(sigma_vals,x); ",
                              "vpc.tw(mean, variance, phi, power); ",
                              "}"
                            ),
                            stop(paste("Gradient for", family, "not implemented."))
  )



  dynamic_func <- eval(parse(text = function_string))
  param_list <- strsplit(param_str, ", ")[[1]]

  # Compute symbolic gradient if method = "symbolic"
  if (method == "symbolic") {
    gradvpc <- sapply(param_list, function(param) {
      deriv_func <- Deriv::Deriv(dynamic_func, param)
      result <- do.call(deriv_func, as.list(param_values))
      result[1]
    })
  }

  # Compute numerical gradient if method = "numerical"
  if (method == "numerical") {
    gradvpc <- numDeriv::grad(function(params) {
      param_list <- as.list(params)
      names(param_list) <- names(param_values)
      do.call(dynamic_func, param_list)
    }, param_values)
    names(gradvpc) <- names(param_values)
  }
  return(gradvpc)
}

gradients <- function(vpcObj, method="numerical") {
  modObj <- vpcObj$modObj
  beta <- modObj$beta
  Sigma <- modObj$Sigma
  phi <- modObj$phi
  family <- modObj$family
  x <- vpcObj$x
  power <- modObj$power
  if (family == "tweedie" && is.null(power)) {
    stop("Parameter 'power' must be provided for Tweedie family in vpcObj.")
  }

  gradient_vpc_engine(beta=beta,
                      Sigma=Sigma,
                      phi=phi,
                      family=family,
                      x=x,
                      power = power,
                      method=method)
}

#' Calculate the Variance-Covariance Matrix for a vpcObj
#'
#' This function computes the variance-covariance matrix for an object of class `vpcObj`.
#'
#' @param object An2 object of class `vpcObj` containing model information and gradients.
#' @param ... Additional arguments (not used in this method).
#'
#' @return The variance of VPC obtained by Delta Method.
#' @export
vcov.vpcObj <- function(object, ...) {
  grad.vpc <- gradients(object)
  vcov.mod <- stats::vcov(object$modObj)
  n <- stats::nobs(object$modObj)
  var.vpc <- (grad.vpc %*% vcov.mod %*% grad.vpc)/n
  return(var.vpc)
}


#' Confidence Intervals for VPC Estimates
#'
#' This function computes confidence intervals for Variance Partition Coefficient (VPC) estimates
#' from an object of class `vpcObj`. The function supports both classical and bootstrap confidence interval types.
#'
#' @param vpcObj An object of class `vpcObj` containing the model fit and VPC estimates.
#' @param alpha Significance level for the confidence intervals. Default is 0.05 (for 95% CI).
#' @param type Character. Specifies the type of confidence interval. Can be either "classical" or "bootstrap".
#'   - "classical": Uses the classical method based on the standard normal distribution.
#'   - "bootstrap": Computes confidence intervals using the parametric bootstrap method.
#' @param iter Integer. The number of bootstrap iterations to perform (for bootstrap type). Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, provides additional information about model convergence and the Hessian matrix's positive definiteness.
#'
#' @return A vector with three elements: Lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   Additional information is included if `verbose = TRUE`. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.vpcObj <- function(vpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap"),
                           iter = 100, num_cores = 1,
                           verbose = FALSE) {
  type <- match.arg(type)

  if (type == "classical") {
    vpc.value <- vpcObj$vpc
    stderr.vpc <- sqrt(stats::vcov(vpcObj))
    crit.val <- stats::qnorm(1 - alpha / 2)
    lwb <- vpc.value - crit.val * stderr.vpc
    upb <- vpc.value + crit.val * stderr.vpc
    result <- c(Lower = lwb, VPC = vpc.value, Upper = upb)

    if (verbose) {
      convergence.code <- vpcObj$modObj$modObj$fit$convergence == 0
      pdHessian <- vpcObj$modObj$modObj$sdr$pdHess
      result <- c(result, list(convergence = convergence.code,
                               pdHess = pdHessian))
    }

    return(result)

  } else if (type == "bootstrap") {
    params <- stats::coef(vpcObj$modObj)
    beta <- params[grep("b", names(params))]
    Sigma_vec <- params[grep("sig", names(params))]
    Sigma <- vec2mat(Sigma_vec)
    family <- vpcObj$modObj$family

    if (family == "negative_binomial") {
      theta <- params[grep("theta", names(params))]
    } else {
      phi <- params[grep("phi", names(params))]
      power <- params[grep("power", names(params))]
    }

    grpVar <- vpcObj$modObj$modObj$modelInfo$grpVar
    frame <- stats::model.frame(vpcObj$modObj)
    ns <- as.numeric(table(frame[, grpVar]))
    X <- frame[, "X"]
    link <- vpcObj$modObj$modObj$modelInfo$family$link
    formula <- formula(vpcObj$modObj$modObj)
    lhs <- as.character(formula[[2]])
    lhs <- gsub("[0-9]+$", "", lhs)
    rhs <- formula[[3]]
    formula <- stats::as.formula(paste(lhs, "~", deparse(rhs)))
    x <- vpcObj$x

    # Simulate new data using the bootstrap method
    data <- glmmVpc::batchGLMMData(beta = beta, ns = ns,
                                   Sigma = Sigma, num = iter,
                                   X = X, theta = theta,
                                   family = family, phi = phi,
                                   power = power, link = link)

    # Fit the model on the simulated data
    fits <- glmmVpc::batchGLMMFit(formula = formula, dataMat = data,
                                  family = family, num_cores = num_cores)

    # Compute the VPC values from the bootstrap fits
    vpcs <- glmmVpc::vpc(fits, x = x)
    vpcs <- sapply(vpcs, function(x) x$vpc)

    # Compute the quantiles to form the confidence interval
    ci <- stats::quantile(vpcs, probs = c(alpha / 2, 1 - alpha / 2))

    result <- c(Lower = ci[1], VPC = stats::median(vpcs), Upper = ci[2])

    return(result)
  }
}


#' Confidence Intervals for VPC Estimates
#'
#' This function computes confidence intervals for Variance Partition Coefficient (VPC) estimates
#' from an object of class `vpcObj`. The function supports both classical and bootstrap confidence interval types.
#'
#' @param VpcObj An object of class `vpcObj` containing the model fit and VPC estimates.
#' @param alpha Significance level for the confidence intervals. Default is 0.05 (for 95% CI).
#' @param type Character. Specifies the type of confidence interval. Can be either "classical" or "bootstrap".
#'   - "classical": Uses the classical method based on the standard normal distribution.
#'   - "bootstrap": Computes confidence intervals using the parametric bootstrap method.
#' @param iter Integer. The number of bootstrap iterations to perform (for bootstrap type). Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, provides additional information about model convergence and the Hessian matrix's positive definiteness.
#'
#' @return A vector with three elements: Lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   Additional information is included if `verbose = TRUE`. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.VpcObj <- function(VpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap"),
                           iter = 100, num_cores = 1,
                           verbose = FALSE) {
  type <- match.arg(type)
  t(sapply(VpcObj, function(vpcObj) stats::confint(vpcObj, alpha=alpha,
                                                   tye=type, iter=iter,
                                                   num_cores=num_cores,
                                                   verbose=verbose)))
}



