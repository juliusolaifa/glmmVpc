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
#' cluster_assignment(c(3, 2, 4))
#'
#' # Generate a factor vector with custom group labels "cluster1", "cluster2"
#' cluster_assignment(c(5, 5), name = "cluster")
cluster_assignment <- function(ns, name="CLUSTER") {
  factor(rep(x=paste0(name, 1:length(ns)), times=ns))
}

#' @export
#' @method as.data.frame glmmDataMatrix
as.data.frame.glmmDataMatrix <- function(x,...) {
  cluster <- colnames(x)
  colnames(x) <- NULL
  x <- unclass(x)
  if (is.matrix(x)) {
    dataMatrix <- t(x)
  } else if (is.vector(x)) {
    dataMatrix <- as.matrix(x)
  } else {
    stop("Unsupported object type for glmmDataMatrix.")
  }
  dataMatrix <- data.frame(dataMatrix)
  dataMatrix$cluster <- cluster
  return(dataMatrix)
}

#' @export
#' @method [ batchglmmDataMatrix
`[.batchglmmDataMatrix` <- function(x, i) {
  if (!inherits(x, "batchglmmDataMatrix")) {
    stop("Object is not of class 'batchglmmDataMatrix'")
  }
  x <- unclass(x)
  x_rows <- x[grepl("^X", rownames(x)), , drop = FALSE]
  if (i > nrow(x) || i < 1) {
    stop("Index out of bounds: no such feature row.")
  }

  feature_row <- x[i+nrow(x_rows), , drop = FALSE]
  result <- rbind(x_rows, feature_row)
  class(result) <- c("glmmDataMatrix", "matrix", "array")
  return(result)
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

#' Extract Predictors (X)
#'
#' This function extracts the predictor variables from a glmmDataMatrix or batchglmmDataMatrix object.
#' If no predictors are found, the function returns NULL.
#'
#' @param object A glmmDataMatrix or batchglmmDataMatrix object.
#'
#' @return A matrix of predictor variables or NULL if no predictors are present.
#' @export
get_x <- function(object) {
  UseMethod("get_x")
}

#' @export
get_x.glmmDataMatrix <- function(object) {
  object <- unclass(object)
  predictors <- grep("^X", rownames(object), value = TRUE)

  if (length(predictors) == 0) {
    return(NULL)
  }
  x <- object[predictors, , drop = FALSE]
  if (nrow(x) == 1) {
    x <- as.vector(x)
  } else {
    x <- t(x)
    rownames(x) <- NULL
  }
  return(x)
}

#' @export
get_x.batchglmmDataMatrix <- get_x.glmmDataMatrix  # Assuming similar structure

#' Extract Response (y)
#'
#' This function extracts the response variable from a glmmDataMatrix or batchglmmDataMatrix object.
#'
#' @param object A glmmDataMatrix or batchglmmDataMatrix object.
#'
#' @return A vector of response variable values.
#' @export
get_y <- function(object) {
  UseMethod("get_y")
}

#' @export
get_y.glmmDataMatrix <- function(object) {
  object <- unclass(object)
  features <- grep("^Feature", rownames(object), value = TRUE)
  y <- object[features, , drop = FALSE]
  return(y)
}

#' @export
get_y.batchglmmDataMatrix <- get_y.glmmDataMatrix  # Assuming similar structure

#' Extract Cluster Information
#'
#' This function extracts the cluster information from the column names of a glmmDataMatrix or batchglmmDataMatrix object.
#'
#' @param object A glmmDataMatrix or batchglmmDataMatrix object.
#'
#' @return A vector of cluster names.
#' @export
get_cluster <- function(object) {
  UseMethod("get_cluster")
}

#' @export
get_cluster.glmmDataMatrix <- function(object) {
  #object <- unclass(object)
  clusters <- as.vector(colnames(object))
  return(clusters)
}

#' @export
get_cluster.batchglmmDataMatrix <- get_cluster.glmmDataMatrix


#' Create a Data Matrix from Covariates, Response Values, and Clusters
#'
#' This function combines a matrix of covariates (X), a matrix of response
#' values (ys),
#' and cluster identifiers into a single data matrix. It ensures that dimensions
#'  are compatible
#' and assigns the appropriate column names.
#'
#' @param X A matrix of covariates. If the number of columns does not match
#' the number of
#'           response values, it will be transposed if the number of rows matches.
#' @param ys A matrix or vector of response values. If a vector, it will be
#' converted to a matrix.
#' @param cluster A character vector containing cluster names to be assigned
#' as column names
#'                for the resulting data matrix.
#'
#' @return A matrix combining covariates, response values, and cluster identifiers.
#'         The resulting matrix has the response values as the last rows and
#'         columns named by the clusters.
#'
#' @export
makeDataMatrix <- function(X, ys, cluster) {
  if (is.null(ys) || is.null(cluster)) {
    stop("Response values (ys) and cluster names must not be NULL.")
  }
  if (!is.null(X) && is.matrix(X)) {
    if (ncol(X) == ncol(ys)) {
      X <- X
    } else if (nrow(X) == ncol(ys)) {
      X <- t(X)
    } else {
      stop("Dimensions of X and ys do not match.")
    }
  }

  result <- rbind(X, ys)
  colnames(result) <- cluster
  class(result) <- c("batchglmmDataMatrix", "matrix", "array")
  attr(result, "num_feat") <- nrow(ys)
  if(is.null(X)) {
    attr(result, "num_covariate") <- 0
  } else if(is.vector(X)) {
    attr(result, "num_covariate") <- 1
  }else{
    attr(result, "num_covariate") <- nrow(X)
  }
  return(result)
}

