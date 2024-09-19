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
