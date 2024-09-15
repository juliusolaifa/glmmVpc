#' Compute the linear predictor and expected response for a GLM/GLMM
#'
#' This function computes the linear predictor (\eqn{\eta = X\beta + ZU}) for a
#' generalized linear model (GLM) or generalized linear mixed model (GLMM), depending
#' on whether the random effects (\eqn{ZU}) are provided. If random effects are not
#' provided, the function defaults to a standard GLM, i.e., \eqn{\eta = X\beta}.
#' It also applies the inverse link function to obtain the expected response \eqn{\mu}.
#'
#' @param X A numeric matrix representing the design matrix for the fixed effects.
#' @param beta A numeric vector representing the fixed effect coefficients.
#' @param Sigma A scalar or symmetric positive-definite matrix representing the
#'              covariance matrix of the distribution.
#' @param ns A vector of group sizes, indicating the number of
#'          observations per group.
#' @param link A character string specifying the link function to use. Supported
#'   options are `"identity"`, `"log"`, `"logit"`, and `"inverse"`.
#' @param warn Logical, controls warning when random effect is not computed.
#'
#' @return A numeric vector representing the expected response \eqn{\mu} after
#'   applying the inverse link function.
#'
#' @examples
#' # Example usage for a GLM with identity link (linear model)
#' X <- matrix(1:10, nrow=10)
#' beta <- c(0.5, 1.5)  # Fixed effect coefficients
#' ns <- c(5,3,2)
#' Sigma <- matrix(c(2,1,1,2),nrow=2)
#' mu <- computeMuGLMM(beta,Sigma, ns, X, link = "identity")
#' print(mu)
#'
#' @export
computeMuGLMM <- function(beta, Sigma=NULL,ns=NULL, X = NULL, link = "identity", warn = TRUE) {
  fixed_effect <- computeFixedEffects(beta=beta, X=X)

  if (is.null(Sigma) || is.null(ns)) {
    random_effect <- 0  # No random effects if either ns or Sigma is NULL
    if (warn) {
      if (is.null(Sigma)) warning("Sigma is NULL, random effects are ignored.")
      if (is.null(ns)) warning("ns is NULL, random effects are ignored.")
    }
  } else {
    random_effect <- generateAndComputeRandomEffects(X = X, ns = ns, Sigma = Sigma)
    if (is.null(random_effect)) random_effect <- 0
  }

  eta <- fixed_effect + random_effect
  mu <- computeMuFromEta(eta, link)
  return(mu)
}

#' Compute mu and sigma for the given design matrices and parameters
#'
#' This helper function computes the linear predictor (mu) and random effect
#' variance (sigma) based on the design matrices for fixed and random effects.
#' If `X` is a vector, it prepends a `1` to represent the intercept and ensures
#' `X` has the same length as `beta`. If `X` is NULL, only the intercept is used.
#' The function handles cases where `Sigma` is either a vector (scalar variance)
#' or a matrix (covariance matrix).
#'
#' @param beta A vector of regression coefficients.
#' @param Sigma A covariance matrix or vector representing variance.
#' @param x  vector for fixed effects. If `x` is NULL, X should only
#' include the intercept.
#'
#' @return A list containing the computed `mu` and `sigma`.
compute_mu_sig <- function(beta, Sigma, x=NULL) {
  if(is.null(beta)) stop("beta can not be NULL")
  if(is.null(Sigma)) stop("Sigma can not be NULL")
  if (is.null(x)) {
    if(length(beta) > 1) {
      stop("beta must be scalar when x is NULL")
    }
    mu <- beta  # Use only the intercept

    if(length(Sigma) > 1) {
      stop("Sigma must be scalar when x is NULL")
    }
    sig <- Sigma

  } else {

      if (!is.vector(x) || !is.numeric(x)) {
          stop("x must be a numeric vector")
      }

      Z <- X <- c(1, x) # Prepend 1 for intercept

      if (length(X) != length(beta)) {
        if(length(beta) == 1) {
          stop(paste("beta is only the intercept", beta,", x should be NULL."))
        }
        stop("Length of X must match the length of beta.")
      }

      mu <- sum(X * beta)

      #Handle univariate and multivariate random effect
      if (is.vector(Sigma) || length(Sigma)==1) {
        sig <- Sigma
      } else if(is.matrix(Sigma)) {
        sig <- as.numeric(Z %*% Sigma %*% Z)
      }else {
        stop("Sigma must be either a scalar (vector of length 1) or a matrix.")
      }
  }

  return(list(mu = mu, sig = sig))
}


#' Custom operator to handle NULL with a default
#'
#' This custom operator returns the left-hand side value unless it is NULL, in
#' which case it returns the right-hand side value.
#'
#' @param a The value to be checked.
#' @param b The default value to return if `a` is NULL.
#'
#' @return Returns `a` if it is not NULL, otherwise returns `b`.
#'
#' @name null_or_default
#' @rdname null_or_default
#'
#' @examples
#' x <- NULL
#' y <- 10
#' x %||% y  # returns 10
#'
#' x <- 5
#' x %||% y  # returns 5
#'
#' @export
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}

