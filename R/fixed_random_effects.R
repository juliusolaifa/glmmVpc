#' Generate Random effects from Multivariate Normal Samples
#'
#' This function generates samples from a multivariate normal distribution with
#' mean 0 and a specified covariance matrix.
#'
#' @param n An integer specifying the number of samples to generate.
#' @param Sigma A scalar or symmetric positive-definite matrix representing the
#'              covariance matrix of the distribution.
#'
#' @return A matrix of dimension \code{n x p}, where \code{p} is the number of
#'          variables (the number of columns in \code{Sigma}). Each row is a
#'          sample from the multivariate normal distribution.
#' @export
#'
#' @examples
#' # Define a covariance matrix
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), ncol=2)
#' # Generate 5 samples from a bivariate normal distribution with 0 vector mean
#' generateRandomEffectsCoefficients(5, Sigma)
generateRandomEffectsCoefficients <- function(n, Sigma) {
  if (is.null(n) || is.null(Sigma)) {
    return(NULL)
  }
  if (!is.numeric(n) || n <= 0 || floor(n) != n) {
    stop("n must be a positive integer.")
  }
  if (is.vector(Sigma) && length(Sigma) == 1) {
    Sigma <- matrix(Sigma, nrow = 1, ncol = 1)
  }
  if (!is.matrix(Sigma) || nrow(Sigma) != ncol(Sigma)) {
    stop("Sigma must be a square matrix.")
  }
  p <- nrow(Sigma)
  tryCatch({
    matrix(MASS::mvrnorm(n=n, mu=rep(0, p), Sigma=Sigma), ncol=p)
  }, error = function(e) {
    stop("Sigma must be a valid positive-definite covariance matrix.")
  })
}

#' Compute Fixed Effects Contribution to the Linear Predictor
#'
#' This function computes the fixed effects part of the
#' linear predictor in a regression model. If the design matrix `X` is provided,
#' it calculates `X %*% beta`. If `X` is `NULL` but `beta` is a scalar
#' (representing an intercept-only model), the function returns a vector
#' of length `n` with the intercept value repeated. The function checks whether
#' an intercept is already included in `X`, and if not, it adds a column of
#' 1's for the intercept.
#'
#' @param X A design matrix (without or with an intercept). If `NULL`,
#' the model is assumed to be intercept-only.
#' @param beta A vector of fixed effect coefficients. If `X` is provided,
#' `beta` should have the same number of elements as columns in `X`. If `X`
#' is `NULL`, `beta` should be a scalar (for an intercept-only model).
#' @param add_intercept Logical to include intercept or not.
#'
#' @return A vector representing the fixed effects contribution to the
#' linear predictor. If `X` is `NULL` and no valid intercept is provided,
#' it returns `0`.
#' @export
#'
#' @examples
#' # Example with X and beta
#' X <- matrix(c(1.5, 2.0, 3.1, 4.0, 5.5, 6.1, 7.3), nrow = 7, ncol = 1)
#' beta <- c(0.5, 1.2)  # Intercept and slope
#' computeFixedEffects(beta, X)
#'
#' # Intercept-only example
#' computeFixedEffects(0.5)
#' computeFixedEffects(beta,X)
computeFixedEffects <- function(beta, X=NULL, add_intercept = TRUE) {
  if (!is.null(beta) && !is.numeric(beta)) {
    stop("beta must be numeric.")
  }

  # For NULL X only beta0 is needed
  if (is.null(X) && !is.null(beta)) {
    if (length(beta) > 1) {
      beta_preview <- paste(utils::capture.output(beta), collapse = "\n")
      stop(sprintf(
        "When X is NULL, beta must have exactly one element (representing the intercept).\nbut beta given as:\n\n%s",
        beta_preview
      ))
    }
    return(beta)
  }

  else if (!is.null(X) && !is.null(beta)) {
    if (!is.numeric(X)) {
      stop("X must be numeric.")
    }
    if(is.vector(X)) {
      X <- matrix(X, ncol=1)
    }
    if (add_intercept && (!all(X[, 1] == 1))) {
      X <- cbind(1, X)
    }

    # Watch out for mismatch, always happen when user do not account for beta0
    if (ncol(X) != length(beta)) {
      X_preview <- paste(utils::capture.output(utils::head(X)), collapse = "\n")
      beta_preview <- paste(utils::capture.output(utils::head(beta)), collapse = "\n")

      stop(sprintf(
        paste(
          "Mismatch between X and beta: X has %d columns, but beta has %d elements.\n",
          "Please ensure that the number of columns in X matches the length
          of beta for matrix multiplication.\n",
          "Here are the first few rows of X and beta:\n\n",
          "X:\n%s\n\n",
          "beta:\n%s"
        ),
        ncol(X), length(beta), X_preview, beta_preview
      ))

    }
    return(X %*% beta)
  } else {
    return(NULL)
  }
}

#' Compute Random Effects Contribution to the Linear Predictor
#'
#' This function computes the random effects part of the linear predictor in a
#' mixed model.
#' It multiplies each design matrix in the list `Z` by the corresponding column
#' in `U` and sums the contributions.
#'
#' @param Z A list of random effects design matrices. Each matrix corresponds
#' to a random effect (e.g., intercept, slope).
#' @param U A matrix of random effect coefficients. Each column corresponds to
#' the coefficients for one random effect.
#'
#' @return A vector representing the random effects contribution to the linear
#' predictor.
#' @export
#'
#' @examples
#' # Example with Z and U
#' ns <- c(2,3,2)
#' X <- rnorm(sum(ns))
#' Z <- generateRandomDesignMatrices(X = X, ns)
#' U <- generateRandomEffectsCoefficients(length(ns), matrix(c(2,1,1,2),2))
#' computeRandomEffects(Z, U)
#' X1 <- matrix(rnorm(sum(ns)*2),ncol=2)
#' Z1 <- generateRandomDesignMatrices(X = X1, ns)
computeRandomEffects <- function(Z, U) {
  if(is.null(Z) && is.null(U)) return(NULL)
  if (!is.list(Z) || !is.matrix(U)) {
    stop("Z must be a list of matrices, and U must be a matrix.")
  }
  # Check if Z is an empty list or U is an empty matrix
  if (length(Z) == 0 || ncol(U) == 0 || nrow(U) == 0) {
    return(NULL)
  }
  if(ncol(U) > length(Z)) {
    stop("Number of random design matrices in Z can not be less than number of columns in U.")
  }
  if(ncol(U) > 1 && length(Z) != ncol(U)) { #ncol(U) = 1 for randpm intercept model
    stop("Number of random design matrices in Z must match number of columns in U in random slope model.")
  }
  random_effect <- rep(0, nrow(Z[[1]]))

  # Loop through each element in Z and the corresponding column in U
  for (i in 1:ncol(U)) { #length(Z) > ncol(U) in g(mu) = alpha +betax + b
      random_effect <- random_effect + Z[[i]] %*% U[, i]
  }

  return(random_effect)
}


#' Generate and Compute Random Effects
#'
#' This function generates the random effects design matrices (`Z`) and random
#' effect coefficients (`U`) based on the input covariate matrix (`X`),
#' sample sizes in each cluster (`ns`), and the variance-covariance matrix
#' of the random effects (`Sigma`). It then computes the random effects contribution
#' to the linear predictor using the `computeRandomEffects` function.
#'
#' @param X A covariate vector or matrix. If it is a vector, it will be used
#' as the covariate for generating random design matrices. If it is a matrix,
#' each column represents a different covariate.
#' @param ns A vector of integers representing the sample sizes in each cluster.
#' @param Sigma A variance-covariance matrix for the random effect coefficients.
#'
#' @return A vector representing the random effects contribution to the linear
#' predictor, computed by multiplying each random effect design matrix (`Z`)
#' with the corresponding column in the random effect coefficient matrix (`U`).
#'
#' @examples
#' # Example usage:
#' X <- rnorm(7)  # Covariate vector
#' ns <- c(2, 3, 2)  # Sample sizes in each cluster
#' Sigma <- diag(2)  # Variance-covariance matrix
#' result <- generateAndComputeRandomEffects(ns, Sigma,X)
#' print(result)
#'
#' @export
generateAndComputeRandomEffects <- function(ns, Sigma, X=NULL) {
  Z <- generateRandomDesignMatrices(ns=ns,X=X)
  U <- generateRandomEffectsCoefficients(n=length(ns), Sigma=Sigma)
  return(computeRandomEffects(Z=Z, U=U))
}
