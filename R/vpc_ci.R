#' Compute Gradient for Variance Partition Coefficient (VPC) Engine
#'
#' This function computes the gradient of the Variance Partition Coefficient (VPC)
#' for different families of distributions, such as Negative Binomial and Tweedie.
#' The gradient can be computed either symbolically or numerically.
#'
#' @param beta A numeric vector of coefficients for the mean function.
#' @param Sigma A numeric matrix representing the covariance structure.
#' @param phi A numeric value representing the dispersion parameter.
#' @param family A character string specifying the family of distributions.
#'   Currently supported families are "negative_binomial" and "tweedie".
#' @param x A numeric vector or matrix of covariates. Default is `NULL`.
#' @param power A numeric value representing the power parameter for the Tweedie family.
#'   Required if `family` is "tweedie". Default is `NULL`.
#' @param method A character string specifying the method to compute the gradient.
#'   Supported methods are "symbolic" and "numerical".
#'
#' @return A numeric vector containing the gradient of the VPC with respect to the parameters.
#'
#' @examples
#' # Example for Negative Binomial family
#' beta <- c(1, 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' phi <- 1.5
#' x <- 1
#' gradient_vpc_engine(beta, Sigma, phi, family = "negative_binomial", x=1,method = "numerical")
#'
#' # Example for Tweedie family
#' power <- 1.5
#' gradient_vpc_engine(beta, Sigma, phi, family = "tweedie", x=1, power = power, method = "numerical")
#'
#' @export
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

  # # Compute symbolic gradient if method = "symbolic"
  # if (method == "symbolic") {
  #   gradvpc <- sapply(param_list, function(param) {
  #     deriv_func <- Deriv::Deriv(dynamic_func, param)
  #     result <- do.call(deriv_func, as.list(param_values))
  #     result[1]
  #   })
  # }

  # Compute numerical gradient if method = "numerical"
  if (method == "numerical") {
    gradvpc <- numDeriv::grad(function(params) {
      param_list <- as.list(params)
      names(param_list) <- names(param_values)
      do.call(dynamic_func, param_list)
    }, param_values)
    names(gradvpc) <- names(param_values)

    hessvpc <- numDeriv::hessian(function(params) {
      param_list <- as.list(params)
      names(param_list) <- names(param_values)
      do.call(dynamic_func, param_list)
    }, param_values)

  }


  return(gradvpc)
}

#' Hessian of the VPC for Negative‑Binomial or Tweedie GLMMs
#'
#' Computes the numerical Hessian of the variance–partitioning coefficient (VPC)
#' with respect to the full parameter vector
#' \eqn{\Theta = (\beta, \phi, \operatorname{vech}(\Sigma))}.
#' The function builds a temporary, family‑specific wrapper around your
#' mean/variance helpers (`mu()`, `sig()`) and VPC calculators
#' (`vpc.nb()`, `vpc.tw()`), then calls **numDeriv**’s
#' \code{\link[numDeriv]{hessian}}.
#'
#' @param beta Numeric vector of fixed‑effect coefficients
#'   \eqn{(\beta_0, \beta_1, \dots)}.
#' @param Sigma Covariance matrix of the random effects (square, symmetric).
#' @param phi  Dispersion parameter:
#'   * negative‑binomial: inverse size (i.e. 1/\eqn{k});
#'   * tweedie: precision \eqn{\phi}.
#' @param family Character string, one of \code{"negative_binomial"} or
#'   \code{"tweedie"}.
#' @param x     Optional design vector (or matrix) of covariates passed to
#'   \code{mu()} and \code{sig()}.  If \code{NULL}, those helpers must handle
#'   the absence of covariates internally.
#' @param power Tweedie power parameter \eqn{p}.  **Required** when
#'   \code{family == "tweedie"}; ignored otherwise.
#' @param method Character, currently only \code{"numerical"} is implemented.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Constructs a parameterised function
#'         \eqn{f(\Theta)} that returns the VPC for the chosen family;
#'   \item Dumps that function into an R expression via \code{eval(parse())};
#'   \item Calls \code{numDeriv::hessian()} at the supplied parameter values.
#' }
#' Helper utilities that must exist in the calling environment:
#' \itemize{
#'   \item \code{construct_beta_param_str()}, \code{construct_sigma_param_str()}
#'         – build comma‑separated parameter names for dynamic evaluation.
#'   \item \code{mat2vec()} – half‑vectorises \eqn{\Sigma}.
#'   \item \code{mu()}, \code{sig()} – mean/variance functions.
#'   \item \code{vpc.nb()}, \code{vpc.tw()} – family‑specific VPC calculators.
#' }
#'
#' @return A square numeric matrix – the Hessian of the VPC with row/column
#'   names matching the parameter vector order
#'   \code{c(beta, phi, vech(Sigma))}.
#'
#' @examples
#' # Example for Negative Binomial family
#' beta <- c(1, 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' phi <- 1.5
#' x <- 1
#'hess_vpc_engine(beta, Sigma, phi, family = "negative_binomial", x=1,method = "numerical")
#'
#' # Example for Tweedie family
#' power <- 1.5
#' hess_vpc_engine(beta, Sigma, phi, family = "tweedie", x=1, power = power, method = "numerical")
#'
#' @export
hess_vpc_engine <- function(beta, Sigma, phi, family, x=NULL, power = NULL, method) {

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

  # Compute numerical gradient if method = "numerical"
  if (method == "numerical") {
    hessvpc <- numDeriv::hessian(function(params) {
      param_list <- as.list(params)
      names(param_list) <- names(param_values)
      do.call(dynamic_func, param_list)
    }, param_values)

    rownames(hessvpc) <- colnames(hessvpc) <- names(param_values)
  }
  return(hessvpc)
}


#' Compute Gradients for a VPC Model
#'
#' This function computes gradients for a given VPC model. It extracts essential
#' parameters from the model object included in \code{vpcObj} and then passes them to the
#' internal function \code{gradient_vpc_engine} for gradient computation.
#'
#' @param vpcObj A list object containing model and data information. It must include:
#'   \describe{
#'     \item{\code{modObj}}{A list with elements \code{beta} (the parameter vector), \code{Sigma}
#'       (covariance matrix), \code{phi} (dispersion parameter), \code{family} (distribution family),
#'       and optionally \code{power} (required when \code{family == "tweedie"}).}
#'     \item{\code{x}}{Predictor variables or design matrix used in the model.}
#'   }
#' @param method A character string specifying the gradient calculation method.
#'   The default is \code{"numerical"}.
#'
#' @return The output of \code{gradient_vpc_engine}, which is the gradient evaluated for the provided
#' model parameters.
#'
#' @details The function checks if the distribution family is \code{"tweedie"}. In that case, it ensures that
#' the \code{power} parameter is provided in \code{modObj}. If \code{power} is \code{NULL}, the function stops
#' and raises an error.
#'
#' @seealso \code{\link{gradient_vpc_engine}}
#'
#' @examples
#' \dontrun{
#'   # Assuming you have a valid vpcObj, you can compute gradients as follows:
#'   result <- gradients(vpcObj, method = "numerical")
#' }
#'
#' @export
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

#' Compute the Hessian for a VPC Model
#'
#' Extracts the model parameters stored in \code{vpcObj} and passes them to
#' the internal helper \code{hess_vpc_engine()} to obtain the numerical (or
#' analytic, when implemented) Hessian of the variance–partitioning
#' coefficient (VPC) with respect to the full parameter vector
#' \eqn{\Theta = (\beta, \phi, \mathrm{vech}(\Sigma))}.
#'
#' @param vpcObj A list containing the fitted model and its data context.
#'   It must include:
#'   \describe{
#'     \item{\code{modObj}}{A sub‑list with elements
#'       \code{beta} (numeric vector of fixed effects),
#'       \code{Sigma} (random‑effect covariance matrix),
#'       \code{phi} (dispersion parameter),
#'       \code{family} (character; currently \code{"negative_binomial"} or
#'       \code{"tweedie"}), and—when \code{family == "tweedie"}—\code{power}.}
#'     \item{\code{x}}{Design matrix or covariate vector used by the model.}
#'   }
#' @param method Character string indicating the differentiation method.
#'   Currently only \code{"numerical"} is available.  Default is
#'   \code{"numerical"}.
#'
#' @return A square numeric matrix: the Hessian of the VPC evaluated at the
#'   supplied parameter values.  Row and column names correspond to
#'   \code{c(beta, phi, vech(Sigma))}.
#'
#' @details
#' For Tweedie models, the power parameter \code{power} must be supplied
#' in \code{modObj}; otherwise the function stops with an error.  All heavy
#' lifting is delegated to \code{hess_vpc_engine()}, which computes the Hessian
#' via \code{numDeriv::hessian()} unless future analytic routines are added.
#'
#' @seealso \code{\link{hess_vpc_engine}}
#'
#' @examples
#' \dontrun{
#'   ## Suppose 'vpcObj' is a previously created VPC fit object
#'   H <- hessians(vpcObj, method = "numerical")
#'   print(H)
#' }
#'
#' @export
hessians <- function(vpcObj, method="numerical") {
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

  hess_vpc_engine(beta=beta,
                      Sigma=Sigma,
                      phi=phi,
                      family=family,
                      x=x,
                      power = power,
                      method=method)
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
    phi <- NULL
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


.make_lower_from_names <- function(nm) {
  lower <- rep(-Inf, length(nm))
  names(lower) <- nm
  # Only sig11 and sig22 are constrained >= 0 (present or absent)
  lower[nm %in% c("sig11", "sig22")] <- 0
  lower
}


#' Sample from a Truncated Multivariate Normal Mixture
#'
#' This function generates random samples from a mixture of truncated multivariate normal
#' distributions. Mixture components are defined by group-wise exclusions of variables.
#' Each component is sampled according to mixture weights (`pis`) with optional truncation
#' to enforce non-negativity constraints on specified variables.
#'
#' @param mean Numeric vector of means for the multivariate normal distribution.
#' @param Sigma Positive-definite covariance matrix, with dimensions matching `mean`.
#' @param pis Numeric vector of mixture probabilities of length \eqn{2^q}, where \eqn{q}
#'   is the number of groups. Must sum to 1 (or will be normalized).
#' @param n Integer, number of samples to generate. Default is 10.
#' @param truncated Logical, whether to truncate variables (default `TRUE`). If `TRUE`,
#'   variables matching `nonneg_prefix` are truncated below at 0.
#' @param groups Named list specifying groups of variable names. Each group defines
#'   coordinates that may be excluded in different mixture components.
#' @param fill_excluded Value used to fill excluded variables (default `0`).
#'
#' @return A numeric matrix of dimension \eqn{n \times p}, where \eqn{p = length(mean)}.
#'   Rows correspond to generated samples, and columns to variables. Excluded coordinates
#'   are filled with `fill_excluded`.
#'
#' @examples
#' set.seed(123)
#' mean  <- c(sig11 = 0, sig22 = 0, x3 = 1)
#' Sigma <- diag(3)
#' pis   <- c(0.5, 0.5)
#' groups <- list(g1 = "sig11")
#'
#' .rtmvnorm_mix_core(mean, Sigma, pis, n = 5, groups = groups)
#'
#' @export
.rtmvnorm_mix_core <- function(mean, Sigma, pis, n = 10, truncated = TRUE,
                                       groups = list(),          # e.g., list(g1="sig11", g2="sig22")
                                       fill_excluded = 0) {

  stopifnot(nrow(Sigma) == length(mean), ncol(Sigma) == length(mean))

  p  <- length(mean)
  nm <- names(mean); if (is.null(nm)) nm <- paste0("x", seq_len(p))
  names(mean) <- nm
  colnames(Sigma) <- rownames(Sigma) <- nm

  # first-version assumption: pis already sum to 1
  stopifnot(abs(sum(pis) - 1) < 1e-10)

  q <- length(groups)
  n_comp <- 2^q
  stopifnot(length(pis) == n_comp)

  # group membership (logical) per group
  group_idx <- lapply(groups, function(g) nm %in% g)

  # mixture allocation (no renormalization)
  w     <- cumsum(pis)
  comps <- findInterval(stats::runif(n), w) + 1L

  # truncation bounds (always truncated per first version)
  lower <- .make_lower_from_names(nm)
  upper <- rep(Inf, p)

  # prefill with zeros for excluded coords
  res <- matrix(fill_excluded, nrow = n, ncol = p, dimnames = list(NULL, nm))

  for (c in seq_len(n_comp)) {
    mask <- comps == c
    if (!any(mask)) next

    # bitmask: bit=1 => EXCLUDE that group's coords (on the boundary face)
    if (q == 0) {
      keep <- rep(TRUE, p)
    } else {
      bits <- as.integer(intToBits(c - 1L))[seq_len(q)] == 1L
      excl <- Reduce(`|`, Map(function(bit, idx) if (bit) idx else FALSE, bits, group_idx), init = FALSE)
      keep <- !excl
    }

    if (!any(keep)) next
    k <- sum(mask)

    # ZERO mean sampling to match first version (do NOT pass mean=)
    res[mask, keep] <- tmvtnorm::rtmvnorm(
      n = k,
      sigma = Sigma[keep, keep, drop = FALSE],
      lower = lower[keep],
      upper = upper[keep]
    )
  }

  res
}


# Public wrappers (clean APIs)

rtmvnorm_mix2q <- function(mean, Sigma, pis, n = 10, truncated = TRUE,
                           exclude1 = "sig11", exclude2 = "sig22",
                           fill_excluded = 0) {
  .rtmvnorm_mix_core(mean, Sigma, pis, n, truncated,
                             groups = list(g1 = exclude1, g2 = exclude2),
                             fill_excluded = fill_excluded)
}

# 2 components (exclude none / exclude one group)
rtmvnorm_mix1q_sig11 <- function(mean, Sigma, pis, n = 10, truncated = TRUE,
                           exclude1 = "sig11",  # default mirrors your branch on sig11
                           fill_excluded = 0) {
  .rtmvnorm_mix_core(mean, Sigma, pis, n, truncated,
                     groups = list(g1 = exclude1),
                     fill_excluded = fill_excluded)
}

# 2 components (exclude none / exclude one group)
rtmvnorm_mix1q_sig22 <- function(mean, Sigma, pis, n = 10, truncated = TRUE,
                           exclude1 = "sig22",  # default mirrors your branch on sig22
                           fill_excluded = 0) {
  .rtmvnorm_mix_core(mean, Sigma, pis, n, truncated,
                              groups = list(g1 = exclude1),
                              fill_excluded = fill_excluded)
}

qmixtnorm.x <- function(mean, Sigma, grad, alpha=0.05, n=1000, thresh=0.01, truncated=TRUE) {
  if(mean["sig11"] < thresh && mean["sig22"] < thresh) {
    mix_func <- rtmvnorm_mix2q
    pis <- c(0.25,0.25,0.25,0.25)
  } else if (mean["sig11"] < thresh) {
    mix_func <- rtmvnorm_mix1q_sig11
    pis <- c(0.5,0.5)
  } else if (mean["sig22"] < thresh) {
    pis <- c(0.5,0.5)
    mix_func <- rtmvnorm_mix1q_sig22
  }
  mixtnorm <- mix_func(mean=mean,Sigma=Sigma,pis=pis,n=n)
  rvpc <- mixtnorm %*% grad
  unname(stats::quantile(rvpc, c(alpha/2,1-alpha/2)))
}

#' Adjusted Confidence Interval Calculation for VPC
#'
#' @param vpcObj An object containing the VPC model with variance-covariance information.
#' @param vpc.value Numeric. The VPC value for which the confidence interval is to be computed.
#' @param forcePD whether the nearest positive matrix should be used
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#' @param n number of bootstrap samples
#' @param thresh Integer. Definition of the boundary
#' @param truncated currently not used
#'
#' @return the confidence interval
#' @export
#'
adjustedc_mixture_ci <- function(vpcObj, vpc.value, forcePD=F, alpha=0.05,
                                 n=1000, thresh=0.01,truncated=TRUE) {
  # vpc.value <- vpcObj$vpc
  fitObj <- vpcObj$modObj
  mean <- stats::coef(fitObj)
  Sigma <- stats::vcov(fitObj, forcePD)
  n.sample <- nobs(fitObj)
  flag <- NULL

  if(mean["sig11"] >= thresh && mean["sig22"] >= thresh) {
    flag <- 0
    cl <- classical_vpc_ci(vpcObj, vpc.value, order=1, alpha = alpha)
    ci <- c(cl[1], cl[2], flag)
    return(ci)
  }

  if(!fitObj$modObj$sdr$pdHess) {
    return(c(NA,NA,NA))
  } else {
    flag <- 1
    grad <- glmmVpc::gradients(vpcObj)
    qmix <- qmixtnorm.x(mean=mean, Sigma=Sigma,
                        grad=grad, alpha=alpha, n=n, thresh=thresh)
    ci <- c(vpc.value - qmix[2], vpc.value - qmix[1], flag)
    return(ci)
  }
}

#' Confidence Interval Calculation for VPC Using Standard Error and Critical Value
#'
#' This function calculates the confidence interval for the Variance Partition Coefficient (VPC)
#' using the standard error and critical value approach. It assumes a normal approximation for the VPC.
#'
#' @param vpcObj An object containing the VPC model with variance-covariance information.
#' @param vpc.value Numeric. The VPC value for which the confidence interval is to be computed.
#' @param forcePD whether the nearest positive matrix should be used
#' @param order An integer indicating if the first order or second order taylo approximation should be used for the delta method
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @export
#'
classical_vpc_ci <- function(vpcObj, vpc.value, forcePD=F, order=1, alpha = 0.05) {
  stderr.vpc <- sqrt(stats::vcov(vpcObj,forcePD, order))
  crit.val <- stats::qnorm(1 - alpha / 2)
  ci <- c(vpc.value - crit.val * stderr.vpc, vpc.value + crit.val * stderr.vpc)

  return(ci)
}



#' Calculate the Variance-Covariance Matrix for a vpcObj
#'
#' This function computes the variance-covariance matrix for an object of class `vpcObj`.
#'
#' @param object An object of class `vpcObj` containing model information and gradients.
#' @param forcePD whether the nearest positive matrix should be used
#' @param order An integer of 1 or 2
#' @param ... Additional arguments (not used in this method).
#'
#' @return The variance of VPC obtained by Delta Method.
#' @export
vcov.vpcObj <- function(object,forcePD=F, order=1, ...) {
  grad.vpc <- gradients(object)
  vcov.mod <- stats::vcov(object$modObj, forcePD)
  var.vpc <- grad.vpc %*% vcov.mod %*% grad.vpc
  if (order == 2) {
    n <- stats::nobs(object$modObj)
    hess.vpc <- hessians(object)
    var.vpc <- var.vpc + (sum(diag(hess.vpc %*% vcov.mod)^2))/(2*n)
  }
  return(var.vpc)
}


#' Confidence Intervals for VPC Estimates
#'
#' This function computes confidence intervals for Variance Partition Coefficient (VPC) estimates
#' from an object of class `vpcObj`. The function supports three confidence interval types: "classical", "bootstrap", and "adjusted".
#'
#' @param vpcObj An object of class `vpcObj` containing the model fit and VPC estimates.
#' @param forcePD whether the nearest positive matrix should be used
#' @param alpha Numeric. The significance level for the confidence intervals. Default is 0.05 (for a 95% CI).
#' @param type Character. Specifies the type of confidence interval to compute. Options include:
#'   - "classical": Based on the standard normal distribution.
#'   - "bootstrap": Using a parametric bootstrap method.
#'   - "adjusted.s": Using the Delta Method, adjusted for boundary conditions (Self & Liang)
#'   - "adjusted.c": Using the Delta Method, adjusted for boundary conditions (Chant)
#' @param thresh Integer. Definition of the boundary
#' @param order An integer indicating if the first order or second order taylo approximation should be used for the delta method
#' @param iter Integer. The number of bootstrap iterations to perform if `type = "bootstrap"`. Default is 100.
#' @param n Integer. The number of monte-carlo samples used for `adjusted_mixture_ci` . Default is 1000.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, additional information is provided regarding model convergence and positive definiteness of the Hessian matrix.
#' @param prob.type Method use in calculating the mixture weights for adjusted method
#'
#' @return A named vector with three elements: the lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   When `verbose = TRUE`, additional diagnostic details on convergence and the Hessian may be included. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.vpcObj <- function(vpcObj, forcePD=F, alpha = 0.05,
                           type = c("classical", "bootstrap",
                                    "adjusted"), thresh=0.01,
                           order=1, num_cores = 4, iter=100,n=1000,
                           verbose = FALSE, prob.type="self") {
  type <- match.arg(type)
  vpc.value <- vpcObj$vpc

  if (type == "classical") {
    ci <- classical_vpc_ci(vpcObj=vpcObj, vpc.value=vpc.value,forcePD=forcePD, order=order,alpha = alpha)
  } else if (type == "bootstrap") {
    ci <- boostrap_vpc_ci(vpcObj, iter = iter, num_cores = num_cores, alpha = alpha)
  } else if (type == "adjusted") {
    ci <- adjustedc_mixture_ci(vpcObj=vpcObj, vpc.value=vpc.value, forcePD=forcePD, alpha = alpha, n = n, thresh=thresh)
  }
  # else if (type == "adjusted.c") {
  #   ci <- adjustedc_mixture_ci(vpcObj, vpc.value, alpha = alpha, n = n,
  #                              truncated = FALSE, prob.type = prob.type)
  # }

  result <- c(Lower = ci[1], VPC = vpc.value, Upper = ci[2])
  names(result) <- c("Lower", "VPC", "Upper")

  if(type == "adjusted") {
    result <- c(Lower = ci[1], VPC = vpc.value, Upper = ci[2],Flag = ci[3])
    names(result) <- c("Lower", "VPC", "Upper","Flag")
  }

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
#' @param forcePD whether the nearest positive matrix should be used
#' @param alpha Significance level for the confidence intervals. Default is 0.05 (for 95% CI).
#' @param type Character. Specifies the type of confidence interval to compute. Options include:
#'   - "classical": Based on the standard normal distribution.
#'   - "bootstrap": Using a parametric bootstrap method.
#'   - "adjusted.s": Using the Delta Method, adjusted for boundary conditions (Self & Liang)
#'   - "adjusted.c": Using the Delta Method, adjusted for boundary conditions (Chant)
#' @param order An integer indicating if the first order or second order taylo approximation should be used for the delta method
#' @param iter Integer. The number of bootstrap iterations to perform (for bootstrap type). Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, provides additional information about model convergence and the Hessian matrix's positive definiteness.
#' @param prob.type Method use in calculating the mixture weights for adjusted method
#' @param thresh Integer. Definition of the boundary
#'
#' @return A vector with three elements: Lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   Additional information is included if `verbose = TRUE`. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.VpcObj <- function(VpcObj, forcePD=F, alpha = 0.05,
                           type = c("classical", "bootstrap",
                                    "adjusted"), thresh=0.01,
                           order=1,iter = 100, num_cores = 4,
                           verbose = FALSE, prob.type = "self") {
  type <- match.arg(type)
  t(sapply(VpcObj, function(vpcObj) stats::confint(vpcObj, forcePD, alpha=alpha,
                                                   type=type, thresh=thresh, order=order,
                                                   iter=iter,num_cores=num_cores,
                                                   verbose=verbose,
                                                   prob.type=prob.type)))
}



# projV <- function(v,idx=NULL) {
#   if(is.null(idx))
#     return(v)
#   else{
#     I <- diag(nrow(v))
#     B <- diag(0,nrow=nrow(v))
#     diag(B)[idx] <- 1
#     P <- I-v%*%B%*%MASS::ginv(t(B)%*%v%*%B)%*%t(B)
#   }
#   V <- P%*%v%*%t(P)
#   return(V[-idx,-idx])
# }
#
#
# adj_ci <- function(vpcObj, N = 1000, seed = NULL) {
#   if (!is.null(seed)) set.seed(seed)
#
#   g <- glmmVpc::gradients(vpcObj)
#   V <- stats::vcov(vpcObj$modObj)
#   idx11 <- which(names(g) == "sig11")
#   idx22 <- which(names(g) == "sig22")
#
#   draw_lin <- function(Vsub, gsub, lower, n) {
#     Z <- tryCatch(tmvtnorm::rtmvnorm(n, sigma = Vsub, lower = lower),
#                   error = function(e) matrix(NA_real_, n, length(gsub)))
#     drop(Z %*% gsub)
#   }
#
#   # ------------------------------------------------------------------
#   draws <- numeric(0)
#
#   if ((idx11) < 0.01 && length(idx22) < 0.01) {
#     # == two‑parameter boundary mixture ==================================
#     V_bi <- V[c(idx11, idx22), c(idx11, idx22)]
#     p1   <- glmmVpc::pr11(vpcObj$modObj)#as.numeric(mvtnorm::pmvnorm(c(0,0), c(Inf,Inf), sigma = V_bi))
#     p2 <- p3 <- 0.25; p4 <- 0.5 - p1
#     stopifnot(p4 >= 0)
#     probs <- c(p1,p2,p3,p4)
#     n_comp <- as.vector(stats::rmultinom(1, N, probs))
#
#     lower1 <- rep(-Inf,length(g)); lower1[c(idx11,idx22)] <- 0
#     if (n_comp[1]) draws <- c(draws, draw_lin(V,g,lower1,n_comp[1]))
#     if (n_comp[2]) {
#       V2 <- projV(V, idx11); g2 <- g[-idx11]
#       lower2 <- rep(-Inf,length(g2)); lower2[names(g2)=="sig22"] <- 0
#       draws <- c(draws, draw_lin(V2,g2,lower2,n_comp[2]))
#     }
#     if (n_comp[3]) {
#       V3 <- projV(V, idx22); g3 <- g[-idx22]
#       lower3 <- rep(-Inf,length(g3)); lower3[names(g3)=="sig11"] <- 0
#       draws <- c(draws, draw_lin(V3,g3,lower3,n_comp[3]))
#     }
#     if (n_comp[4]) {
#       V4 <- projV(V, c(idx11,idx22)); g4 <- g[-c(idx11,idx22)]
#       draws <- c(draws, draw_lin(V4,g4,rep(-Inf,length(g4)),n_comp[4]))
#     }
#   } else {
#     # == single‑parameter boundary mixture ================================
#     idx <- if (idx11 < 0.01) idx11 else idx22
#     par_name <- names(g)[idx]
#     n_comp <- as.vector(stats::rmultinom(1, N, prob = c(0.5,0.5)))
#     lower_full <- rep(-Inf,length(g)); lower_full[idx] <- 0
#     if (n_comp[1]) draws <- c(draws, draw_lin(V,g,lower_full,n_comp[1]))
#     if (n_comp[2]) {
#       V2 <- projV(V, idx); g2 <- g[-idx]
#       lower2 <- rep(-Inf,length(g2)); lower2[names(g2)==par_name] <- 0
#       draws <- c(draws, draw_lin(V2,g2,lower2,n_comp[2]))
#     }
#   }
#
#   delta <- stats::quantile(draws, c(0.975,0.025), na.rm = TRUE)
#   vpc   <- vpcObj$vpc
#   c(lower = vpc - delta[1], upper = vpc - delta[2])
# }
#
