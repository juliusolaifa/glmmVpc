#' Extract model parameters by family
#'
#' This function extracts the relevant model parameters (e.g., beta, Sigma, phi, power)
#' based on the specified family from a fitted GLMM model object.
#'
#' @param family A character string representing the family (e.g., "nb" for Negative Binomial).
#' @param modObj The fitted GLMM model object.
#'
#' @return A list of extracted model parameters, for random intercept we return
#' the variance and not standard deviation.
#' @export
extractParametersByFamily <- function(family, modObj) {

  # Helper function to extract common parameters
  extract_common_params <- function(modObj) {
    beta <- unname(glmmTMB::fixef(modObj)$cond)
    Sigma <- unname(glmmTMB::VarCorr(modObj)$cond)[[1]]
    # Strip all attributes from Sigma
    attr(Sigma, "stddev") <- attr(Sigma, "correlation") <- attr(Sigma, "blockCode") <- NULL
    # If Sigma is scalar, square it to represent the variance
    if(length(Sigma) == 1) {
      Sigma <- Sigma^2
    }

    return(list(beta = beta, Sigma = Sigma))
  }

  # Switch case to handle specific families
  switch(family,
         "negative_binomial" = {
           params <- extract_common_params(modObj)
           phi <- glmmTMB::sigma(modObj)  # theta equivalent
           c(params, list(phi = phi, family = family))
         },
         "tweedie" = {
           params <- extract_common_params(modObj)
           phi <- glmmTMB::sigma(modObj)
           power <- unname(glmmTMB::family_params(modObj))
           c(params, list(phi = phi, power = power, family = family))
         },
         "gaussian" = {
           params <- extract_common_params(modObj)
           sigma_e <- glmmTMB::sigma(modObj)
           c(params, list(sigma_e=sigma_e, family = family))
         },
         stop(paste("Unsupported family:", family))
  )
}



#' Fit a single GLMM model
#'
#' This function fits a GLMM using the specified formula, family, and data,
#' and extracts the relevant model parameters in a way that mirrors data generation by family.
#'
#' @param formula A formula describing the model.
#' @param data A data frame containing the variables in the model.
#' @param family A character string representing the family (e.g., "nb" for Negative Binomial).
#'
#' @return A list with model parameters.
#' @export
singleGLMMFit <- function(formula, data, family) {

  if(!is.data.frame(data)) {
    data <- as.data.frame(data)
  }

  if (missing(formula) || missing(data) || missing(family)) {
    stop("Formula, data, and family are required inputs.")
  }
  glmmTMBfamily <- tryCatch({
    get_glmmTMBfamily(family)
  }, error = function(e) {
    stop("Error in creating glmmTMB family: ", e$message)
  })

  modObj <- tryCatch({
    glmmTMB::glmmTMB(formula = formula, data = data, family = glmmTMBfamily)
  }, error = function(e) {
    message("Error fitting model: ", e$message)
    return(NULL)
  })
  if (is.null(modObj)) {
    return(NULL)
  }
  params <- extractParametersByFamily(family, modObj)
  sample_size <- nrow(data)
  var_cov_matrix <- stats::vcov(modObj, full=TRUE)

  params$sample_size <- sample_size
  params$var_cov_matrix <- var_cov_matrix
  params$modObj <- modObj


  class(params) <- "glmmfit"
  return(params)
}

nobs.glmmfit <- function(object, ...) {
  if (!is.null(object$sample_size)) {
    return(object$sample_size)
  } else {
    stop("Sample size (nobs) is not available in the glmmfit object.")
  }
}

vcov.glmmfit <- function(object, ...) {
  if (!is.null(object$var_cov_matrix)) {
    return(object$var_cov_matrix)
  } else {
    stop("Variance-covariance matrix (vcov) is not available in the glmmfit object.")
  }
}

model.frame.glmmfit <- function(object, ...) {
  if (!is.null(object$modObj)) {
    return(stats::model.frame(object$modObj, ...) && inherits(object$modObj, "glmmTMB"))
  } else {
    stop("Model object (modObj) is not available in the glmmfit object.")
  }
}

logLik.glmmfit <- function(object, ...) {
  if (!is.null(object$modObj) && inherits(object$modObj, "glmmTMB")) {
    return(stats::logLik(object$modObj, ...))
  } else {
    stop("Model object (modObj) is not available or not a valid glmmTMB object in the glmmfit object.")
  }
}


#' Batch Generalized Linear Mixed Model (GLMM) Fit for Multiple Features
#'
#' This function fits a GLMM for each feature in the provided data matrix.
#' The formula for each feature is dynamically generated.
#'
#' @param formula A formula specifying the model to fit, such as
#' Feature ~ 1 + (1|cluster).
#' @param dataMat A batchglmmDataMatrix containing both covariates and
#' feature data. Covariates should match `cov.pattern`.
#' @param family A description of the error distribution and link function to
#' be used in the model (e.g., "gaussian", "poisson").
#' @param cov.pattern A regular expression pattern to identify covariates in
#' the `dataMat`. Default is "^X".
#' @param num_cores Number of computer cores to use in training.
#'
#' @return A list of model fits, one for each feature.
#' @export
batchGLMMFit <- function(formula, dataMat, family, cov.pattern = "^X", num_cores = 1) {

  max_cores <- parallel::detectCores() - 1

  if (num_cores < 1) {
    warning("num_cores must be at least 1. Setting num_cores to 1 (no parallelism).")
    num_cores <- 1
  } else if (num_cores > max_cores) {
    warning(paste("num_cores exceeds the available cores. Using", max_cores, "cores instead."))
    num_cores <- max_cores
  }

  print(paste("Fitting GLMM for ", family))
  num_feat <- attr(dataMat, "num_feat")
  fits <- pbapply::pblapply(1:num_feat, function(i) {
    data <- as.data.frame(dataMat[i])
    formula_string <- paste(paste0(deparse(formula[[2]]),i), "~", deparse(formula[[3]]))
    dynamic_formula <- stats::as.formula(formula_string)
    fit <- tryCatch({
      singleGLMMFit(formula = dynamic_formula, data = data, family = family)
    }, error = function(e) {
      warning(paste("Error fitting model for Feature",i, ":", e$message))
      return(NULL)
    })

    return(fit)
  }, cl = num_cores)

  fits <- Filter(Negate(is.null), fits)
  return(fits)
}


#' Convert family name to glmmTMB family object
#'
#' This function converts a string representing a GLMM family (e.g., "nb" for Negative Binomial)
#' into a corresponding glmmTMB family object.
#'
#' @param family A string representing the family (e.g., "nb", "poisson", "gaussian").
#' @return A glmmTMB family object.
#' @export
get_glmmTMBfamily <- function(family) {
  switch(family,
         "gaussian" = stats::gaussian,
         "negative_binomial" = glmmTMB::nbinom2,
         "tweedie" = glmmTMB::tweedie,
         stop("Unsupported family")
  )
}


