#' Generate Data Based on Specified Distribution Family
#'
#' This function generates random data from a specified distribution family.
#' Supported families include Gaussian, Poisson, Negative Binomial, Tweedie,
#' Gamma, and Binomial.
#'
#' @param family A string specifying the distribution family. Supported values
#' are: "gaussian", "poisson", "negative_binomial", "tweedie", "gamma",
#' "binomial".
#' @param mu A numeric value representing the mean or rate
#' (depending on the family).
#' @param ns A vector of integers representing the number of samples
#' for each group.
#' @param args A list of additional arguments specific to each distribution
#' family. The required
#' arguments for each family are:
#'   - "gaussian": `sd` (standard deviation, optional, defaults to 1)
#'   - "negative_binomial": `theta` (shape parameter, required)
#'   - "tweedie": `p` (power parameter, required) and `phi`
#'   (dispersion parameter, required)
#'   - "gamma": `shape` (shape parameter, required)
#'
#' @return A vector of randomly generated data based on the specified family
#' and parameters.
#'
#' @examples
#' # Generate Gaussian data with mean 0 and standard deviation 2
#' generateDataByFamily("gaussian", mu = 0, ns = c(10, 10, 10),
#' args = list(sd = 2))
#'
#' # Generate Poisson data with rate 5
#' generateDataByFamily("poisson", mu = 5, ns = c(10, 10, 10), args = list())
#'
#' # Generate Negative Binomial data
#' generateDataByFamily("negative_binomial", mu = 10, ns = c(10, 10),
#' args = list(theta = 1.5))
#'
#' @export
generateDataByFamily <- function(family, mu, ns, args) {
  if (!is.numeric(mu)) stop("'mu' must be a numeric vector.")
  if (!is.numeric(ns) || any(ns <= 0)) stop("'ns' must be a vector of positive integers.")
  if (!is.list(args)) stop("'args' must be a list.")

  switch(family,
         "gaussian" = stats::rnorm(sum(ns), mean = mu, sd = args$sd %||% 1),
         "poisson" = stats::rpois(sum(ns), lambda = mu),
         "negative_binomial" = mapply(MASS::rnegbin, n = 1, mu = mu, theta = args$theta),
         "tweedie" = mapply(tweedie::rtweedie, n=1, mu = mu, power = args$power, phi = args$phi),
         "gamma" = stats::rgamma(sum(ns), shape = args$shape, rate = args$shape / mu),
         "binomial" = stats::rbinom(sum(ns), size = 1, prob = mu),
         stop("Unsupported family")
  )
}

#' Generate GLMM Data for a Single Model
#'
#' This function generates data for a single generalized linear mixed model (GLMM).
#' It computes the mean of the model, generates response data using the specified
#' distribution family, and returns the data along with the design matrix.
#'
#' @param X A matrix of predictors (design matrix) for the fixed effects.
#' @param beta A vector of fixed effect coefficients.
#' @param Sigma A covariance matrix for the random effects.
#' @param ns A vector of sample sizes for each group.
#' @param family A string specifying the distribution family (default is "gaussian").
#'               Supported families include "gaussian", "poisson", "negative_binomial", "tweedie",
#'               "gamma", and "binomial".
#' @param link A string specifying the link function (default is "identity").
#' @param ... Additional arguments to be passed to the family-specific data generation functions.
#'
#' @return A matrix where rows represent the predictors and the response variable (named "Feature").
#'
#' @examples
#' X <- matrix(1:30)
#' beta <- c(0.5, 0.2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' ns <- c(10, 10, 10)
#' singleGLMMData(beta, ns, Sigma, X, family = "gaussian", link = "identity")
#'
#' @export
singleGLMMData <- function(beta, ns, Sigma=NULL, X=NULL, family = "gaussian",
                           link = "identity", ...) {

  args <- list(...)

  family_args <- validate_family_params(family, args)
  if (is.null(X) && !is.null(Sigma) && (length(Sigma) != 1) ){
      stop(paste("For random intercept-only models (X = NULL), Sigma must be a scalar.\n",
                 "Current Sigma is a matrix:\n",
                 paste(utils::capture.output(Sigma), collapse = "\n")))
  }
  mu <- computeMuGLMM(X=X, beta=beta, Sigma=Sigma, ns=ns, link = link)
  y <- tryCatch({
    generateDataByFamily(family, mu, ns, family_args)
  }, error = function(e) {
    warning("Error in generating data for family '", family, "': ",
            conditionMessage(e))
    return(NULL)
  })

  if (is.null(y)) {
    warning("Data generation failed, returning NULL")
    return(NULL)
  }

  # If X is not NULL, handle the design matrix
  if (!is.null(X)) {
    if(is.vector(X)) X <- matrix(X, ncol=1)
    Xt <- if (all(X[, 1] == 1) && ncol(X) > 1) t(as.matrix(X[, -1])) else t(X)
    if(nrow(Xt)==1) rownames(Xt) <-  "X" else  rownames(Xt) <-  paste0("X", 1:nrow(Xt))
    y <- matrix(unlist(y), nrow = 1)
    rownames(y) <- "Feature"
    data <- rbind(Xt, y)
  } else {
    # If X is NULL, return only the generated response y
    data <- matrix(unlist(y), nrow=1)
    rownames(data) <- "Feature"
  }
  colnames(data) <- cluster_assignment(ns)
  class(data) <- c("glmmDataMatrix", "matrix", "array")
  return(data)
}

#' Generate GLMM Data for Multiple Features
#'
#' This function generates data for multiple features using generalized linear
#' mixed models (GLMMs).
#' It calls the `singleGLMMData` function for each feature and combines the
#' results into a matrix.
#'
#' @param num The number of features (Samples) to generate.
#' @param X A matrix of predictors (design matrix) for the fixed effects.
#' @param beta A vector of fixed effect coefficients.
#' @param Sigma A covariance matrix for the random effects.
#' @param ns A vector of sample sizes for each group.
#' @param family A string specifying the distribution family
#' (default is "gaussian").Supported families include "gaussian", "poisson",
#' "negative_binomial", "tweedie", "gamma", and "binomial".
#' @param link A string specifying the link function (default is "identity").
#' @param ... Additional arguments to be passed to the family-specific data
#' generation functions.
#'
#' @return A matrix where rows represent the predictors (design matrix `X`)
#' and generated features. The feature rows are named as "Feature 1",
#' "Feature 2", etc.
#'
#' @examples
#' X <- matrix(1:10, nrow = 10)
#' beta <- c(0.5, 0.2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' ns <- c(5,3,2)
#' batchGLMMData(beta = beta, ns = ns, Sigma = Sigma, num = 5, X = X,
#'               family = "gaussian", link = "identity")
#'
#' @export
batchGLMMData <- function(beta, ns, Sigma=NULL, num=1, X=NULL,
                          family = "gaussian", link = "identity",...) {
  if (num <= 0) {
    warning("num_feat must be positive. Returning NULL.")
    return(NULL)
  }
  dataMatrix <- lapply(1:num, function(x) {
    single_data <- singleGLMMData(beta=beta, ns=ns,Sigma=Sigma, X=X, family=family, link=link, ...)
    if (!is.null(single_data)) {
      return(single_data[nrow(single_data), , drop = FALSE])
    } else {
      warning("Feature generation failed for iteration ", x)
      return(NULL)
    }
  })

  dataMatrix <- Filter(Negate(is.null), dataMatrix)

  if (length(dataMatrix) == 0) {
    warning("No data was successfully generated, returning NULL")
    return(NULL)
  }

  dataMatrix <- do.call(rbind, dataMatrix)

  if (!is.null(X)) {
    Xt <- t(X)
    if (nrow(Xt) == 1) {
      rownames(Xt) <- "X"
    } else {
      rownames(Xt) <- paste0("X", 1:nrow(Xt))
    }
    rownames(dataMatrix) <- paste0("Feature", 1:nrow(dataMatrix))
    fullData <- rbind(Xt, dataMatrix)
  } else {
    # Handle the case where X is NULL (intercept-only model)
    rownames(dataMatrix) <- paste0("Feature", 1:nrow(dataMatrix))
    fullData <- dataMatrix
  }

  class(fullData) <- c("batchglmmDataMatrix", "matrix", "array")
  attr(fullData, "num_feat") <- num
  attr(fullData, "num_covariate") <- ifelse(is.null(X), 0, nrow(Xt))
  return(fullData)
}


validate_family_params <- function(family, args) {
  required_params <- list(
    "negative_binomial" = c("theta"),
    "tweedie" = c("power", "phi"),
    "gamma" = c("shape"),
    "poisson" = c(),
    "gaussian" = c(),
    "binomial" = c()
  )

  # Check if the family is supported
  if (!family %in% names(required_params)) {
    stop(paste("Unsupported family:", family))
  }


  required <- required_params[[family]]
  # Find missing parameters
  missing_params <- setdiff(required, names(args))
  if (length(missing_params) > 0) {
    stop(paste("Missing required parameters for", family,
               "family:", paste(missing_params, collapse = ", ")))
  }


  # Extract family-specific parameters from args
  family_specific_args <- as.list(args[required])
  return(family_specific_args)
}
