gradient_vpc_engine <- function(beta, Sigma, phi, family, x=NULL, p = NULL, method="symbolic") {

  if (family == "tweedie" && is.null(p)) {
    stop("Parameter 'p' is required for Tweedie family.")
  }

  # Construct the function string dynamically
  beta_str <- construct_beta_param_str(beta)
  Sigma_str <- construct_sigma_param_str(Sigma)
  param_str <- switch(family,
                      "negative_binomial" = paste(beta_str, "phi", Sigma_str, sep = ", "),
                      "tweedie" = paste(beta_str, "phi", "p", Sigma_str, sep = ", ")
  )

  param_values <- switch(family,
                         "negative_binomial" = c(
                           stats::setNames(beta, paste0("b", seq_along(beta) - 1)),
                           phi = phi,
                           mat2vec(Sigma)),
                         "tweedie" = c(
                           stats::setNames(beta, paste0("b", seq_along(beta) - 1)),
                           phi = phi,
                           p = p,
                           mat2vec(Sigma))
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
                              "vpc.tw(mean, variance, phi, p); ",
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

gradients <- function(vpcObj, method="symbolic") {
  beta <- attr(vpcObj, "beta")
  Sigma <- attr(vpcObj, "Sigma")
  phi <- attr(vpcObj, "phi")
  family <- attr(vpcObj, "family")
  x <- attr(vpcObj, "x")
  p <- attr(vpcObj, "p")
  if (family == "tweedie" && is.null(p)) {
    stop("Parameter 'p' must be provided for Tweedie family in vpcObj.")
  }

  gradient_vpc_engine(beta=beta,
                      Sigma=Sigma,
                      phi=phi,
                      family=family,
                      x=x,
                      p = p,
                      method=method)
}

confint.vpc <- function(vpcObj, alpha=0.05) {
  vcov_glmm <- attr(vpcObj, "vcov")
  sample_size <- attr(vpcObj, "n")
  crit_val <- stats::qnorm(1-alpha/2)

  grad_vpc <- gradients(vpcObj)

  variance_vpc <- (grad_vpc %*% vcov_glmm %*% grad_vpc)/sample_size #Delta Theoerem

  stderr_vpc <- sqrt(variance_vpc)

  vpc_value <- vpcObj$vpc_value

  # Compute confidence intervals
  lower_bound <- vpc_value - crit_val * stderr_vpc
  upper_bound <- vpc_value + crit_val * stderr_vpc

  # Return a named vector with the VPC value and confidence intervals
  return(c(VPC = vpc_value, Lower = lower_bound, Upper = upper_bound))

}
