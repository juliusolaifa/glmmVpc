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

#' Generate random samples from a mixture of (optionally truncated) multivariate normals
#'
#' This function generates \code{n} random draws from a mixture of up to four different
#' (truncated or non-truncated) multivariate normal components, determined by
#' \code{pis}. The function internally checks for entries named \code{"sig11"} or
#' \code{"sig22"} in \code{mean} to decide how to truncate each component.
#'
#' @param mean A numeric vector of means for each dimension. The names of \code{mean} are used
#'   to determine which dimensions potentially get truncated (\code{"sig11"} or \code{"sig22"}).
#' @param Sigma A covariance matrix for the underlying multivariate normal distribution.
#' @param pis A numeric vector of mixture proportions. Must sum up to \code{1}, and must have
#'   length up to \code{4} to match the four possible components in the mixture.
#' @param n An integer specifying the number of random draws to generate.
#' @param truncated A logical flag. If \code{TRUE}, samples for certain components
#'   are drawn from a truncated multivariate normal using \pkg{tmvtnorm}; if \code{FALSE},
#'   the function uses \pkg{mvtnorm} without truncation for those components.
#'
#' @details
#' The function does the following for each draw:
#' \enumerate{
#'   \item Draws a uniform random number \code{u}.
#'   \item Checks which mixture component \code{u} falls into (based on \code{pis}).
#'   \item Depending on the component, it selects which subset of dimensions to truncate at 0
#'         (if \code{truncated=TRUE}) or sample normally:
#'         \itemize{
#'           \item Component 1: No truncation; draw from \code{mvtnorm::rmvnorm(1, Sigma)}.
#'           \item Component 2: Truncate dimensions \emph{not} labeled \code{"sig11"}.
#'           \item Component 3: Truncate dimensions \emph{not} labeled \code{"sig22"}.
#'           \item Component 4: Neither \code{"sig11"} nor \code{"sig22"} dimensions are truncated.
#'         }
#' }
#'
#' @return A matrix with \code{n} rows and \code{length(mean)} columns, where each row is
#'   one random draw from the specified mixture distribution.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#'
#' # Example: Suppose we have two dimensions, named "sig11" and "otherDim".
#' # We'll draw from a mixture with two components (50/50) where
#' #   - First component: standard MVN(0, I)
#' #   - Second component: truncated in dimension "otherDim" >= 0
#'
#' my_mean  <- c(sig11 = 0, otherDim = 0)
#' my_Sigma <- diag(2)
#' my_pis   <- c(0.5, 0.5)  # two-component mixture
#'
#' samples <- rmixtnorm(mean = my_mean, Sigma = my_Sigma, pis = my_pis,
#'                      n = 1000, truncated = TRUE)
#' head(samples)
#' }
#'
#' @export
rmixtnorm <- function(mean, Sigma, pis, n=10, truncated=TRUE) {
  # Generate n uniform samples for selecting mixture components
  rng <- stats::runif(n)
  pis <- cumsum(pis)

  # Identify truncation conditions
  has_sig11 <- names(mean) %in% "sig11"
  has_sig22 <- names(mean) %in% "sig22"
  has_sig11_sig22 <- names(mean) %in% c("sig11", "sig22")
  lower <- c(-Inf, -Inf, 0, 0, -Inf, 0)
  # Determine which component each sample belongs to
  idxs <- findInterval(rng, pis) + 1  # Assigns component indices (1,2,3,4)

  # Pre-allocate result matrix
  result <- matrix(0, nrow=n, ncol=length(mean))
  colnames(result) <- names(mean)

  # Sample from different components **in batches** (instead of looping)

  ## Component 1: No truncation, full MVN sample
  mask1 <- (idxs == 1)
  if(any(mask1)) {
    lower1 <- lower
    result[mask1,] <- tmvtnorm::rtmvnorm(sum(mask1), sigma=Sigma, lower=lower1)
  }

  ## Component 2: Truncate dimensions **not labeled** "sig11"
  mask2 <- (idxs == 2)
  if(any(mask2)) {
    idx <- !has_sig11
    lower2 <- lower[idx]
    result[mask2, idx] <- tmvtnorm::rtmvnorm(sum(mask2),
                                             sigma=Sigma[idx, idx, drop=FALSE],
                                             lower=lower2)
  }

  ## Component 3: Truncate dimensions **not labeled** "sig22"
  mask3 <- (idxs == 3)
  if(any(mask3)){
    idx <- !has_sig22
    lower3 <- lower[idx]
    result[mask3, idx] <- tmvtnorm::rtmvnorm(sum(mask3),
                                             sigma=Sigma[idx, idx, drop=FALSE],
                                             lower=lower3)
  }
  ## Component 4: No truncation, subset of MVN
  mask4 <- (idxs == 4)
  if(any(mask4)){
    idx <- !has_sig11_sig22
    lower4 <- lower[idx]
    result[mask4, idx] <- tmvtnorm::rtmvnorm(sum(mask4),
                                             sigma=Sigma[idx, idx, drop=FALSE],
                                             lower=lower4)
  }

  return(result)
}

# rmixtnorm <- function(mean, Sigma, pis, n=10, truncated=TRUE) {
#   # Generate n uniform samples for selecting mixture components
#   rng <- stats::runif(n)
#   pis <- cumsum(pis)
#
#   # Identify truncation conditions
#   has_sig11 <- names(mean) %in% "sig11"
#   has_sig22 <- names(mean) %in% "sig22"
#   has_sig11_sig22 <- names(mean) %in% c("sig11", "sig22")
#
#   # Determine which component each sample belongs to
#   idxs <- findInterval(rng, pis) + 1  # Assigns component indices (1,2,3,4)
#
#   # Pre-allocate result matrix
#   result <- matrix(0, nrow=n, ncol=length(mean))
#   colnames(result) <- names(mean)
#
#   # Sample from different components **in batches** (instead of looping)
#
#   ## Component 1: No truncation, full MVN sample
#   mask1 <- (idxs == 1)
#   result[mask1, ] <- mvtnorm::rmvnorm(sum(mask1), sigma=Sigma)
#
#   ## Component 2: Truncate dimensions **not labeled** "sig11"
#   mask2 <- (idxs == 2)
#   if (any(mask2)) {
#     idx <- !has_sig11
#     lower <- ifelse(has_sig22[idx], 0, -Inf)
#     result[mask2, idx] <- if (truncated) {
#       tmvtnorm::rtmvnorm(sum(mask2), sigma=Sigma[idx, idx, drop=FALSE], lower=lower)
#     } else {
#       mvtnorm::rmvnorm(sum(mask2), sigma=Sigma[idx, idx, drop=FALSE])
#     }
#   }
#
#   ## Component 3: Truncate dimensions **not labeled** "sig22"
#   mask3 <- (idxs == 3)
#   if (any(mask3)) {
#     idx <- !has_sig22
#     lower <- ifelse(has_sig11[idx], 0, -Inf)
#     result[mask3, idx] <- if (truncated) {
#       tmvtnorm::rtmvnorm(sum(mask3), sigma=Sigma[idx, idx, drop=FALSE], lower=lower)
#     } else {
#       mvtnorm::rmvnorm(sum(mask3), sigma=Sigma[idx, idx, drop=FALSE])
#     }
#   }
#
#   ## Component 4: No truncation, subset of MVN
#   mask4 <- (idxs == 4)
#   if (any(mask4)) {
#     idx <- !has_sig11_sig22
#     result[mask4, idx] <- mvtnorm::rmvnorm(sum(mask4), sigma=Sigma[idx, idx, drop=FALSE])
#   }
#
#   return(result)
# }

# rmixtnorm <- function(mean, Sigma, pis, n=10, truncated = TRUE) {
#   rng <- stats::runif(n)
#   pis <- cumsum(pis)
#   result <- matrix(NA,nrow=n, ncol = length(mean))
#   colnames(result) <- names(mean)
#   has_sig11 <- names(mean) %in% "sig11"
#   has_sig22 <- names(mean) %in% "sig22"
#   has_sig11_sig22 <- names(mean) %in% c("sig11","sig22")
#
#   sample_component <- function(sigma, lower, use_trunc) {
#     if (use_trunc) {
#       return(tmvtnorm::rtmvnorm(1, sigma = sigma, lower = lower))
#     } else {
#       return(mvtnorm::rmvnorm(1, sigma = sigma))
#     }
#   }
#
#   for(i in 1:n) {
#     dat <- numeric(length(mean))
#     names(dat) <- names(mean)
#     if(rng[i] <= pis[1]) {
#       dat <- mvtnorm::rmvnorm(1,sigma = Sigma)
#     }
#     else if(rng[i] <= pis[2]) {
#       idx <- !has_sig11
#       lower = ifelse(has_sig22[idx], 0, -Inf)
#       dat.temp <- sample_component(Sigma[idx,idx, drop=F],lower, truncated)
#       dat[idx] <- dat.temp
#     }
#     else if(rng[i] <= pis[3]) {
#       idx <- !has_sig22
#       lower = ifelse(has_sig11[idx], 0, -Inf)
#       dat.temp <- sample_component(Sigma[idx,idx, drop=F],lower, truncated)
#       dat[!has_sig22] <- dat.temp
#     }
#     else {
#       idx <- !has_sig11_sig22
#       dat.temp <- mvtnorm::rmvnorm(1,sigma=Sigma[idx, idx])
#       dat[idx] <- dat.temp
#     }
#     result[i,] <- dat
#   }
#   return(result)
# }

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
#' @param truncated Logical. Whether truncation should be used.
#'
#' @return A numeric vector of length 2 containing the lower and upper quantiles for the confidence interval, calculated as \code{c(alpha/2, 1 - alpha/2)}.
#'
#' @export
#'
qmixtnorm <- function(mean, Sigma, pis, grad, alpha=0.05, n=1000, truncated=TRUE) {
  mixtnorm <- rmixtnorm(mean=mean, Sigma=Sigma, pis=pis, n=n, truncated=truncated)
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

#' Confidence Interval Computation for VPC using Mixture Normal Quantiles
#'
#' This function calculates the confidence interval for the Variance Partition Coefficient (VPC)
#' by computing quantiles from a mixture of normal distributions. The function extracts
#' necessary parameters from a fitted model object and applies the `qmixtnorm` function.
#'
#' @param vpcObj An object containing the fitted model and additional VPC-related information.
#' @param vpc.value Numeric. The VPC value for which the confidence interval is to be computed.
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#' @param n Integer. Number of samples to draw for quantile estimation (default is 1000).
#' @param truncated Logical. Whether truncation should be used.
#' @param prob.type Method use in calculating the mixture weights
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @export
#'
adjustedc_mixture_ci <- function(vpcObj, vpc.value, alpha=0.05, n=1000,
                                 truncated=TRUE,
                                 prob.type="self") {
  # Extract fitted model from vpcObj
  fitObj <- vpcObj$modObj

  # Obtain mean vector and covariance matrix from the model
  mean <- stats::coef(fitObj)
  Sigma <- stats::vcov(fitObj)
  n.sample <- nobs(fitObj)

  # Calculate mixture proportions
  if(!fitObj$modObj$sdr$pdHess) {
    return(c(NA,NA))
  }else{
    p11 <- pr11(fitObj, prob.type)
    pis <- c(p11, 0.25, 0.25, 0.5 - p11)

    # Get gradient vector for transformation
    grad <- gradients(vpcObj)

    # Compute the confidence interval using mixture normal quantiles
    qmix <- qmixtnorm(mean=mean, Sigma=Sigma, pis=pis,
                      grad=grad, alpha=alpha, n=n, truncated=truncated)
    ci <- c(vpc.value - qmix[2], vpc.value - qmix[1])
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
#' @param order An integer indicating if the first order or second order taylo approximation should be used for the delta method
#' @param alpha Numeric. The significance level for the confidence interval (default is 0.05).
#'
#' @return A numeric vector of length 2 containing the lower and upper bounds of the confidence interval.
#' @export
#'
classical_vpc_ci <- function(vpcObj, vpc.value, order=1, alpha = 0.05) {
  stderr.vpc <- sqrt(stats::vcov(vpcObj,order))
  crit.val <- stats::qnorm(1 - alpha / 2)
  ci <- c(vpc.value - crit.val * stderr.vpc, vpc.value + crit.val * stderr.vpc)

  return(ci)
}



#' Calculate the Variance-Covariance Matrix for a vpcObj
#'
#' This function computes the variance-covariance matrix for an object of class `vpcObj`.
#'
#' @param object An object of class `vpcObj` containing model information and gradients.
#' @param order An integer of 1 or 2
#' @param ... Additional arguments (not used in this method).
#'
#' @return The variance of VPC obtained by Delta Method.
#' @export
vcov.vpcObj <- function(object,order=1, ...) {
  grad.vpc <- gradients(object)
  vcov.mod <- stats::vcov(object$modObj)
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
#' @param alpha Numeric. The significance level for the confidence intervals. Default is 0.05 (for a 95% CI).
#' @param type Character. Specifies the type of confidence interval to compute. Options include:
#'   - "classical": Based on the standard normal distribution.
#'   - "bootstrap": Using a parametric bootstrap method.
#'   - "adjusted.s": Using the Delta Method, adjusted for boundary conditions (Self & Liang)
#'   - "adjusted.c": Using the Delta Method, adjusted for boundary conditions (Chant)
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
confint.vpcObj <- function(vpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap",
                                    "adjusted.s", "adjusted.c"),
                           num_cores = 4, iter=100,n=1000,
                           verbose = FALSE, prob.type="self") {
  type <- match.arg(type)
  vpc.value <- vpcObj$vpc

  if (type == "classical") {
    ci <- classical_vpc_ci(vpcObj, vpc.value, alpha = alpha)
  } else if (type == "bootstrap") {
    ci <- boostrap_vpc_ci(vpcObj, iter = iter, num_cores = num_cores, alpha = alpha)
  } else if (type == "adjusted.s") {
    ci <- adjustedc_mixture_ci(vpcObj, vpc.value, alpha = alpha, n = n)
  } else if (type == "adjusted.c") {
    ci <- adjustedc_mixture_ci(vpcObj, vpc.value, alpha = alpha, n = n,
                               truncated = FALSE, prob.type = prob.type)
  }

  result <- c(Lower = ci[1], VPC = vpc.value, Upper = ci[2])
  names(result) <- c("Lower", "VPC", "Upper")

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
#'   - "adjusted.s": Using the Delta Method, adjusted for boundary conditions (Self & Liang)
#'   - "adjusted.c": Using the Delta Method, adjusted for boundary conditions (Chant)
#' @param iter Integer. The number of bootstrap iterations to perform (for bootstrap type). Default is 100.
#' @param num_cores Integer. The number of cores to use for parallel computation in the bootstrap method. Default is 1.
#' @param verbose Logical. If `TRUE`, provides additional information about model convergence and the Hessian matrix's positive definiteness.
#' @param prob.type Method use in calculating the mixture weights for adjusted method
#'
#' @return A vector with three elements: Lower bound, the VPC estimate, and the upper bound of the confidence interval.
#'   Additional information is included if `verbose = TRUE`. Returns `NA` if the model did not converge or if the Hessian is not positive definite.
#'
#' @export
confint.VpcObj <- function(VpcObj, alpha = 0.05,
                           type = c("classical", "bootstrap",
                                    "adjusted.s", "adjusted.c"),
                           iter = 100, num_cores = 4,
                           verbose = FALSE, prob.type = "self") {
  type <- match.arg(type)
  t(sapply(VpcObj, function(vpcObj) stats::confint(vpcObj, alpha=alpha,
                                                   type=type, iter=iter,
                                                   num_cores=num_cores,
                                                   verbose=verbose,
                                                   prob.type=prob.type)))
}



projV <- function(v,idx=NULL) {
  if(is.null(idx))
    return(v)
  else{
    I <- diag(nrow(v))
    B <- diag(0,nrow=nrow(v))
    diag(B)[idx] <- 1
    P <- I-v%*%B%*%MASS::ginv(t(B)%*%v%*%B)%*%t(B)
  }
  V <- P%*%v%*%t(P)
  return(V[-idx,-idx])
}



