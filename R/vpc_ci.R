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

rmixtnorm <- function(mean, Sigma, pis, n=10) {
  rng <- stats::runif(n)
  print(pis)
  pis <- cumsum(pis)
  result <- matrix(NA,nrow=n, ncol = length(mean))
  colnames(result) <- names(mean)
  has_sig11 <- names(mean) %in% "sig11"
  has_sig22 <- names(mean) %in% "sig22"
  has_sig11_sig22 <- names(mean) %in% c("sig11","sig22")
  for(i in 1:n) {
    dat <- numeric(length(mean))
    names(dat) <- names(mean)
    if(rng[i] <= pis[1]) {
      dat <- mvtnorm::rmvnorm(1,mean=mean,sigma = Sigma)
    }
    else if(rng[i] > pis[1] && rng[i] <= pis[2]) {
      lower = ifelse(has_sig22[!has_sig11], 0, -Inf)
      truncnorm_moments <- tmvtnorm::mtmvnorm(mean = mean[!has_sig11],
                                              sigma = Sigma[!has_sig11,!has_sig11],
                                              lower = lower)

      dat.temp <- tmvtnorm::rtmvnorm(1,mean=truncnorm_moments$tmean,
                                     sigma=truncnorm_moments$tvar,lower=lower)
      dat[!has_sig11] <- dat.temp
    }
    else if(rng[i] > pis[2] && rng[i] <= pis[3]) {
      lower = ifelse(has_sig11[!has_sig22], 0, -Inf)
      truncnorm_moments <- tmvtnorm::mtmvnorm(mean = mean[!has_sig22],
                                              sigma = Sigma[!has_sig22,!has_sig22],
                                              lower = lower)
      dat.temp <- tmvtnorm::rtmvnorm(1,mean=truncnorm_moments$tmean,
                                     sigma=truncnorm_moments$tvar,lower=lower)
      dat[!has_sig22] <- dat.temp
    }
    else {
      lower <- ifelse(has_sig11_sig22[!has_sig11_sig22], 0, -Inf)
      dat.temp <- tmvtnorm::rtmvnorm(1,mean=mean[!has_sig11_sig22],
                                     sigma=Sigma[!has_sig11_sig22, !has_sig11_sig22],
                                     lower=lower)
      dat[!has_sig11_sig22] <- dat.temp
    }
    result[i,] <- dat
  }
  result
  return(result)
}

#' Quantile Computation for a Mixture of Multivariate & Truncated Normals
#'
#' This function calculates the quantiles for a mixture of multivariate normal
#' distributions
#' based on the provided mean vector, covariance matrix, mixture proportions,
#' and gradient vector.
#' The quantiles are returned for the confidence interval specified by `alpha`.
#'
#' @param mean A numeric vector representing the mean of each component in the multivariate normal distribution.
#' @param Sigma A covariance matrix associated with the multivariate normal distribution.
#' @param pis A numeric vector of mixing proportions for the components of the mixture model. The elements of `pis` should sum to 1.
#' @param grad A numeric vector representing the gradient, used to transform the sampled mixture data.
#' @param alpha A numeric value specifying the significance level for the confidence interval (default is 0.05, which corresponds to a 95% confidence interval).
#' @param n An integer specifying the number of samples to draw from the mixture distribution (default is 100).
#'
#' @return A numeric vector of length 2 containing the lower and upper quantiles for the confidence interval, calculated as \code{c(alpha/2, 1 - alpha/2)}.
#'
#' @export
#'
qmixtnorm <- function(mean, Sigma, pis, grad, alpha=0.05, n=100) {
  mixtnorm <- rmixtnorm(mean=mean, Sigma=Sigma, pis=pis, n=n)
  rvpc <- mixtnorm %*% grad
  stats::quantile(rvpc, c(alpha/2, 1-alpha/2))
}

#' Bootstrap Confidence Interval for VPC
#'
#' This function calculates a bootstrap confidence interval for the Variance Partition Coefficient (VPC) by
#' simulating new data, fitting the model, and computing the VPC values from the bootstrap fits.
#'
#' @param vpcObj An object containing model information, including coefficients, family, and grouping variable.
#' @param iter Integer. Number of iterations for the bootstrap simulation.
#' @param num_cores Integer. Number of cores to use for parallel processing during model fitting.
#' @param alpha Numeric. Significance level for the confidence interval (default is 0.05).
#'
#' @return A numeric vector of length 2 containing the lower and upper quantiles for the confidence interval.
#' @export
#'
boostrap_vpc_ci <- function(vpcObj, iter = 100, num_cores = 4, alpha = 0.05) {
  # Extract model parameters
  params <- stats::coef(vpcObj$modObj)
  beta <- params[grep("b", names(params))]
  Sigma_vec <- params[grep("sig", names(params))]
  Sigma <- vec2mat(Sigma_vec)
  family <- vpcObj$modObj$family

  # Check for family-specific parameters
  if (family == "negative_binomial") {
    theta <- params[grep("theta", names(params))]
  } else {
    phi <- params[grep("phi", names(params))]
    power <- params[grep("power", names(params))]
  }

  # Extract additional model info
  grpVar <- vpcObj$modObj$modObj$modelInfo$grpVar
  frame <- stats::model.frame(vpcObj$modObj)
  ns <- as.numeric(table(frame[, grpVar]))
  X <- frame[, "X"]
  link <- vpcObj$modObj$modObj$modelInfo$family$link

  # Construct formula for model fitting
  formula <- formula(vpcObj$modObj$modObj)
  lhs <- as.character(formula[[2]])
  lhs <- gsub("[0-9]+$", "", lhs)
  rhs <- formula[[3]]
  formula <- stats::as.formula(paste(lhs, "~", deparse(rhs)))
  x <- vpcObj$x

  # Simulate new data using bootstrap method
  data <- glmmVpc::batchGLMMData(beta = beta, ns = ns,
                                 Sigma = Sigma, num = iter,
                                 X = X, theta = theta,
                                 family = family, phi = phi,
                                 power = power, link = link)

  # Fit the model on simulated data
  fits <- glmmVpc::batchGLMMFit(formula = formula, dataMat = data,
                                family = family, num_cores = num_cores)

  # Compute VPC values from bootstrap fits
  vpcs <- glmmVpc::vpc(fits, x = x)
  vpcs <- sapply(vpcs, function(x) x$vpc)

  # Compute confidence interval
  ci <- stats::quantile(vpcs, probs = c(alpha / 2, 1 - alpha / 2))

  return(ci)
}

#' Confidence Interval Computation for VPC using Mixture Normal Quantiles
#'
#' This function calculates the confidence interval for the Variance Partition Coefficient (VPC)
#' by computing quantiles from a mixture of normal distributions. The function extracts
#' necessary parameters from a fitted model object and applies the `qmixtnorm` function.
#'
#' @param vpcObj An object containing the fitted model and additional VPC-related information.
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#' @param n Integer. Number of samples to draw for quantile estimation (default is 1000).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @export
#'
adjustedc_mixture_ci <- function(vpcObj, alpha = 0.05, n = 1000) {
  # Extract fitted model from vpcObj
  fitObj <- vpcObj$modObj

  # Obtain mean vector and covariance matrix from the model
  mean <- stats::coef(fitObj)
  Sigma <- stats::vcov(fitObj)

  # Calculate mixture proportions
  p11 <- pr11(fitObj)
  pis <- c(p11, 0.25, 0.25, 0.5 - p11)

  # Get gradient vector for transformation
  grad <- gradients(vpcObj)

  # Compute the confidence interval using mixture normal quantiles
  ci <- qmixtnorm(mean = mean, Sigma = Sigma, pis = pis, grad = grad, alpha = alpha, n = n)

  return(ci)
}

#' Confidence Interval Calculation for VPC Using Standard Error and Critical Value
#'
#' This function calculates the confidence interval for the Variance Partition Coefficient (VPC)
#' using the standard error and critical value approach. It assumes a normal approximation for the VPC.
#'
#' @param vpcObj An object containing the VPC model with variance-covariance information.
#' @param vpc.value Numeric. The VPC value for which the confidence interval is to be computed.
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @export
#'
classical_vpc_ci <- function(vpcObj, vpc.value, alpha = 0.05) {
  # Calculate standard error of the VPC
  stderr.vpc <- sqrt(stats::vcov(vpcObj))

  # Determine critical value based on significance level
  crit.val <- stats::qnorm(1 - alpha / 2)

  # Compute the confidence interval
  ci <- c(vpc.value - crit.val * stderr.vpc, vpc.value + crit.val * stderr.vpc)

  return(ci)
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
#' from an object of class `vpcObj`. The function supports three confidence interval types: "classical", "bootstrap", and "adjusted".
#'
#' @param vpcObj An object of class `vpcObj` containing the model fit and VPC estimates.
#' @param alpha Numeric. The significance level for the confidence intervals. Default is 0.05 (for a 95% CI).
#' @param type Character. Specifies the type of confidence interval to compute. Options include:
#'   - "classical": Based on the standard normal distribution.
#'   - "bootstrap": Using a parametric bootstrap method.
#'   - "adjusted": Using the Delta Method, adjusted for boundary conditions.
#' @param iter Integer. The number of bootstrap iterations to perform if `type = "bootstrap"`. Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, additional information is provided regarding model convergence and positive definiteness of the Hessian matrix.
#'
#' @return A named vector with three elements: the lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   When `verbose = TRUE`, additional diagnostic details on convergence and the Hessian may be included. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.vpcObj <- function(vpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap", "adjusted"),
                           num_cores = 1, iter=100,verbose = FALSE) {
  type <- match.arg(type)
  vpc.value <- vpcObj$vpc

  if (type == "classical") {
    ci <- classical_vpc_ci(vpcObj, vpc.value, alpha = 0.05)
  } else if (type == "bootstrap") {
    ci <- boostrap_vpc_ci(vpcObj, iter = iter, num_cores = 4, alpha = alpha)
  } else if (type == "adjusted") {
    ci <- adjustedc_mixture_ci(vpcObj, alpha = alpha, n = 1000)
  }

  result <- c(Lower = ci[1], VPC = vpc.value, Upper = ci[2])

  if (verbose) {
    convergence.code <- vpcObj$modObj$modObj$fit$convergence == 0
    pdHessian <- vpcObj$modObj$modObj$sdr$pdHess
    result <- c(result, list(convergence = convergence.code,
                             pdHess = pdHessian))
  }
  return(result)
}


#' Confidence Intervals for VPC Estimates
#'
#' This function computes confidence intervals for Variance Partition Coefficient (VPC) estimates
#' from an object of class `vpcObj`. The function supports both classical and bootstrap confidence interval types.
#'
#' @param VpcObj An object of class `vpcObj` containing the model fit and VPC estimates.
#' @param alpha Significance level for the confidence intervals. Default is 0.05 (for 95% CI).
#' @param type Character. Specifies the type of confidence interval to compute. Options include:
#'   - "classical": Based on the standard normal distribution.
#'   - "bootstrap": Using a parametric bootstrap method.
#'   - "adjusted": Using the Delta Method, adjusted for boundary conditions.
#' @param iter Integer. The number of bootstrap iterations to perform (for bootstrap type). Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, provides additional information about model convergence and the Hessian matrix's positive definiteness.
#'
#' @return A vector with three elements: Lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   Additional information is included if `verbose = TRUE`. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.VpcObj <- function(VpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap", "adjusted"),
                           iter = 100, num_cores = 1,
                           verbose = FALSE) {
  type <- match.arg(type)
  t(sapply(VpcObj, function(vpcObj) stats::confint(vpcObj, alpha=alpha,
                                                   type=type, iter=iter,
                                                   num_cores=num_cores,
                                                   verbose=verbose)))
}



