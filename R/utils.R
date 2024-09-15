#' Linear Predictor (Eta)
#'
#' This function calculates the linear predictor (often denoted as  \eqn{\eta}) in
#' a generalized linear model (GLM) or linear mixed model (LMM). The linear
#' predictor is computed as a combination of fixed effects  \eqn{\beta} and
#' random effects  \eqn{Z}.
#'
#' @param X A numeric matrix of fixed effect covariates. Each row represents an
#' observation, and each column represents a covariate.
#' @param beta A numeric vector of fixed effect coefficients. The length of
#' `beta` should match the number of columns in `X`.
#' @param Z A numeric matrix of random effect covariates. Each row represents
#' an observation, and each column represents a random effect covariate.
#' @param u A numeric vector of random effect coefficients. The length of `u`
#' should match the number of columns in `Z`.
#'
#' @return A numeric vector representing the linear predictor, with each element
#' corresponding to an observation.
#'
#' @export
#'
#' @examples
#' X <- matrix(c(1, 1, 1, 0, 1, 2), nrow=3, ncol=2)
#' beta <- c(0.5, 2)
#' Z <- matrix(c(1, 0, 1, 1, 1, 0), nrow=3, ncol=2)
#' u <- c(1.5, -0.5)
#' eta(X, beta, Z, u)
eta <- function(X, beta, Z, u) {
  X %*% beta + Z %*% u
}


#' Generate Group Factor Labels
#'
#' This function generates a factor vector with group labels, where each group
#' is labeled with a specified name and a unique identifier
#' (e.g., "group1", "group2", etc.). The number of times each group label
#' appears is determined by the corresponding value in the `ns` vector.
#'
#' @param ns A numeric vector where each element specifies the number of
#' repetitions for the corresponding group.
#' @param name A character string specifying the base name for the group labels.
#' The default is "group".
#'
#' @return A factor vector with group labels repeated according to the values in
# `ns`.
#' @export
#'
#' @examples
#' # Generate a factor vector with group labels "group1", "group2", "group3"
#' groups(c(3, 2, 4))
#'
#' # Generate a factor vector with custom group labels "cluster1", "cluster2"
#' groups(c(5, 5), name = "cluster")
groups <- function(ns, name="group") {
  factor(rep(x=paste0(name, 1:length(ns)), times=ns))
}


cluster_assignment <- function(ns, name="CLUSTER") {
  factor(rep(x=paste0(name, 1:length(ns)), times=ns))
}

as_dataframe <- function(dataMatrix) {
  if(is.matrix(dataMatrix)) {
    dataMatrix <- t(dataMatrix)
  } else if(is.vector(dataMatrix)){
    dataMatrix <- as.matrix(dataMatrix)
  }
  dataframe <- as.data.frame(dataMatrix)
  dataframe$cluster <- rownames(dataMatrix)
  rownames(dataframe) <- NULL
  return(dataframe)
}


#' Generate Random Binary Vectors with Equal Probability
#'
#' This function generates random binary vectors (composed of 0s and 1s) for
#' each element in the input vector `ns`. Each 0 and 1 is equally likely, with
#' the length of each binary vector determined by the corresponding element in
#' `ns`.
#'
#' @param ns A numeric vector where each element specifies the length of the
#' binary vector to be generated.
#'
#' @return A concatenated vector of 0s and 1s, where the length of each section
#' of the vector corresponds to the respective value in `ns`.
#' @export
#'
#' @examples
#' # Generate binary vectors of lengths 3, 5, and 2
#' rgen01(c(3, 5, 2))
#'
#' # Generate binary vectors of lengths 10, 20, and 30
#' rgen01(c(10, 20, 30))
rgen01 <- function(ns) {
  unlist(
    lapply(ns, function(x)
      sample(rep_len(c(0, 1), length.out = x))
    )
  )
}


#' Convert Vector to Symmetric Matrix
#'
#' This function converts a vector into a symmetric matrix by filling the lower
#' triangular part of the matrix with the elements of the vector and mirroring
#' it to the upper triangular part. The size of the matrix is automatically
#' determined based on the length of the vector.
#'
#' @param vec A numeric vector containing the elements to be placed in the lower
#' triangular part of the matrix. The length of `vec` should be compatible with
#' forming a symmetric matrix.
#'
#' @return A symmetric numeric matrix with dimensions determined by the length
#' of `vec`.
#'
#' @export
#'
#' @examples
#' # Convert a vector to a symmetric 3x3 matrix
#' vec <- c(2,1,2)
#' vec2mat(vec)
#'
#' # Example with a longer vector for a larger matrix
#' vec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' vec2mat(vec)
vec2mat <- function(vec) {
  k <- length(vec)
  n <- (-1 + sqrt(1 + 8 * k)) / 2
  mat <- matrix(0, nrow = n, ncol = n)
  mat[!upper.tri(mat)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}

mat2vec <- function(mat) {
  n <- dim(mat)[1]
  mat <- mat[!upper.tri(mat)]
  names(mat) <- unlist(sapply(1:n, function(k) {
    paste0('sig',k,k:n)
  }))
  return(mat)
}

#' GLMMTMB Familiies
#'
#' This function returns the appropriate family function for `glmmTMB` based on the input string.
#'
#' @param family A string indicating the family to be used. Available options are:
#'               "nb" (negative binomial), "tw" (Tweedie), "ga" (Gaussian).
#'
#' @return The corresponding family function for `glmmTMB`.
#' @export
#'
#' @examples
#' glmmTMBfamily("nb")
#' glmmTMBfamily("tw")
#' glmmTMBfamily("ga")
glmmTMBfamily <- function(family){
    switch(family,
            nb=glmmTMB::nbinom2,
            tw=glmmTMB::tweedie,
            ga=stats::gaussian,
            stop(family, "is not a valid family.
                     Available options are: nb, tw, ga.")
    )
}


vcov.fit <- function(mod) {
  vcov_ <- stats::vcov(mod, full=TRUE)
  rownames(vcov_) <- colnames(vcov_) <- attr(rownames(vcov_), "names")
  vcov_
}


#' Compute the k-th Moment of a Log-Normal Distribution
#'
#' This function computes the k-th moment of a log-normal distribution given the mean `mu`,
#' variance `Sigma` (or sigma^2), and a constant `k`.
#'
#' @param mu Numeric value representing the mean of the underlying normal distribution.
#' @param Sigma Numeric value representing the variance (sigma^2) of the underlying normal distribution.
#' @param k Numeric value representing the moment to compute.
#'
#' @return Returns the k-th moment of the log-normal distribution.
#' @export
#'
#' @examples
#' logNormMoment(0, 1, 1) # First moment (mean) of a log-normal distribution
#' logNormMoment(0.5, 2, 2) # Second moment of a log-normal distribution
logNormMoment <- function(mu, Sigma, k) {
  return(exp(k*mu + k^2*Sigma/2))
}


#' Compute moments of the log-normal distribution
#'
#' This function computes the moments (e.g., first, second, power-th moments)
#' of a log-normal distribution
#' based on the mean (`mu`) and variance (`sigma^2`) of the corresponding
#' normal distribution.
#'
#' @param mu The mean (linear predictor) on the log scale.
#' @param variance The variance (random effect variance) on the log scale.
#' @param power The power for higher-order moments (e.g., for Tweedie models).
#'
#' @return A list containing the first moment (`Ex`), the second moment
#' (`Ex2`), and, if applicable, the power-th moment (`Ex_p`).
compute_lnmoments <- function(mu, variance, power = 1) {
  Ex <- exp(mu + 0.5 * variance)
  Ex2 <- exp(2 * mu + 2 * variance)
  Ex_p <- if (power != 1) exp(power * mu + 0.5 * power^2 * variance) else Ex

  return(list(Ex = Ex, Ex2 = Ex2, Ex_p = Ex_p))
}



