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


mvn <- function(n, Sigma) {
  p <- dim(Sigma)[1]
  MASS::mvrnorm(n=n, mu=rep(0,p), Sigma=Sigma)
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
#' This function converts a vector into a symmetric matrix by filling the lower triangular part of the matrix with the elements of the vector and mirroring it to the upper triangular part. The size of the matrix is automatically determined based on the length of the vector.
#'
#' @param vec A numeric vector containing the elements to be placed in the lower triangular part of the matrix. The length of `vec` should be compatible with forming a symmetric matrix.
#'
#' @return A symmetric numeric matrix with dimensions determined by the length of `vec`.
#' @export
#'
#' @examples
#' # Convert a vector to a symmetric 3x3 matrix
#' vec <- c(1, 2, 3, 4, 5, 6)
#' vec2mat(vec)
#'
#' # Example with a longer vector for a larger matrix
#' vec <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' vec2mat(vec)
vec2mat <- function(vec) {
  n <- as.integer(stats::uniroot(function(x) {x^2 + x - 2*length(vec)},
                                 lower = 1, upper = 10)$root)
  mat <- matrix(0, nrow = n, ncol = n)
  mat[!upper.tri(mat)] <- vec
  mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]
  mat
}

