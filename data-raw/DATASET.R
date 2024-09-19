# Code to prepare `paramNBmat` dataset goes here

# Parameters for Negative Binomial
b0 <- c(3, 5)
b1 <- c(5, 7)
sig11 <- 2
sig12 <- 1
sig22 <- 2
theta <- seq(0.01, 0.4, 0.02)

# Create a parameter matrix
paramNBmat <- expand.grid(b0 = b0, b1 = b1,
                          sig11 = sig11, sig12 = sig12,
                          sig22 = sig22, theta = theta)

#' Parameter Matrix for Negative Binomial
#'
#' A dataset containing combinations of parameters for a generalized linear mixed model.
#'
#' @format A parameter matrix for Negative Binomial for generating glmmDataMatrix:
#' \describe{
#'   \item{b0}{Values for parameter b0.}
#'   \item{b1}{Values for parameter b1.}
#'   \item{sig11}{Variance for the random intercept.}
#'   \item{sig12}{Covariance between the random intercept and random slope.}
#'   \item{sig22}{Variance for the random slope.}
#'   \item{theta}{Dispersion parameter for the Negative Binomial model.}
#' }
#' @source Generated for use in GLMM fitting.
"paramNBmat"

# Save the dataset
usethis::use_data(paramNBmat, overwrite = TRUE)

# Example for generating 'X' (ensure this function is defined)
ns <- c(5, 3, 2)
X <- rgen01(ns)  # Make sure this function is defined in your package

# Generate batch data using the parameter matrix and X
glmmDataMatrixNB <- batchGLMMDataUsingMatrix(paramMat = paramNBmat, ns = ns, X = X,
                                       family = "negative_binomial", link = "log")

#' Batch Data for Negative Binomial GLMM
#'
#' A dataset generated for a generalized linear mixed model using Negative Binomial parameters.
#'
#' @format A matrix representing the generated batch data:
#' \describe{
#'   \item{X}{Covariate matrix.}
#'   \item{Feature}{Generated response variable based on the model.}
#' }
#' @source Simulated data for demonstration purposes.
"glmmDataMatrixNB"

# Save the batch data
usethis::use_data(glmmDataMatrixNB, overwrite = TRUE)
