#' Generate Random Intercept Design Matrix
#'
#' This function generates a design matrix for random intercepts, where each
#' group has its own random intercept, and all observations within a group
#' share the same random effect.
#'
#' @param ns A vector specifying the number of observations per group
#'        (cluster sizes).
#'
#' @return A matrix where each column represents a group, and each row
#'          corresponds to an observation. The entries are 1 if the
#'          observation belongs to the corresponding group and 0 otherwise.
#' @export
#'
#' @examples
#' ns <- c(3, 3)  # Two groups, each with 3 observations
#' Z <- generateRandomInterceptMatrix(ns)
#' print(Z)
generateRandomInterceptMatrix <- function(ns) {
  if(is.null(ns)) return (NULL)
  if (!is.numeric(ns) || any(ns <= 0)) {
    stop("ns must be a numeric vector of positive integers.")
  }
  num_clusters <- length(ns)
  total_obs <- sum(ns)

  if (num_clusters == 0) return(NULL)

  Z0 <- matrix(0, nrow = total_obs, ncol = num_clusters)
  group_indices <- rep.int(1:num_clusters, ns)
  Z0[cbind(seq_len(total_obs), group_indices)] <- 1
  return(Z0)
}

#' Generate Random Design Matrices for Mixed Models
#'
#' This function generates random design matrices for
#' random intercepts and random slopes. If `X` contains a column of 1's as the
#' first column (indicating an intercept), the intercept matrix is not generated
#' separately. If `X` is NULL, only the intercept matrix is generated.
#'
#' @param X Optional. A matrix or vector of fixed effect predictors.
#' The first column is assumed to be an intercept if all values are 1.
#' @param ns A vector of group sizes, indicating the number of
#' observations per group.
#'
#' @return A list containing the intercept and slope matrices. If `X` is NULL,
#' only the intercept matrix is returned.
#'
#' @export
#'
#' @examples
#' ns <- c(2, 3, 2)
#' X <- matrix(c(rep(1,7),seq(2,20,length.out=7)), nrow = 7, ncol = 2)
#' print(X)
#' generateRandomDesignMatrices(ns,X)
generateRandomDesignMatrices <- function(ns,X = NULL) {
  if (is.null(ns)) return(NULL)
  if (!is.numeric(ns) || any(ns <= 0)) {
    stop("ns must be a numeric vector of positive integers.")
  }
  Z0 <- generateRandomInterceptMatrix(ns)
  Z_list <- list(Z0)
  if (!is.null(X) && length(X) > 0) {

    if (is.vector(X)) {
      X <- matrix(X, ncol = 1)
    }
    if (nrow(X) != sum(ns)) {
      stop("The number of rows in X must match the total number of observations in ns.")
    }

    if (all(X[, 1] == 1)) {
      X <- X[, -1, drop = FALSE]
    }
    if(ncol(X) >= 1) {
      for (j in 1:ncol(X)) {
        Z_list[[length(Z_list) + 1]] <- Z0 * X[, j]
      }
    }
  }
  return(Z_list)
}
