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

  if(modObj$fit$convergence == 1 || modObj$sdr$pdHess) {
    print("Re-fitting")
    modObj <- stats::update(modObj, control=glmmTMB::glmmTMBControl(
                                                  optimizer=stats::optim,
                                                  optArgs=list(method="BFGS")))
  }

  if (is.null(modObj)) {
    return(NULL)
  }
  params <- extractParametersByFamily(family, modObj)
  params$modObj <- modObj


  class(params) <- "glmmfit"
  return(params)
}

#' @importFrom stats nobs
#' @export
#' @method nobs glmmfit
nobs.glmmfit <- function(object, ...) {
    return(stats::nobs(object$modObj))
}


#' Variance-Covariance Matrix for GLMM Fit
#'
#' This function computes the Asymptotic variance-covariance matrix of the
#' regression parameters, family specific parameters and variance component of the random effect
#'
#' @param object A fitted GLMM object of class \code{glmmfit}.
#' @param ... Additional arguments (currently not used).
#'
#' @return A matrix containing the variance-covariance matrix of the estimated parameters.
#' The row and column names are set to the parameter names based on the family of the model.
#'
#' @examples
#' \dontrun{
#'   fit <- glmmfit(data, family = "binomial")
#'   vcov_matrix <- vcov(fit)
#' }
#'
#' @seealso \code{\link[stats]{vcov}} for the generic variance-covariance method.
#'
#' @export
vcov.glmmfit <- function(object, ...) {
  modObj <- object$modObj
  vcovObj <- stats::vcov(modObj, full = TRUE)
  rownames(vcovObj) <- colnames(vcovObj) <- par_names(object, object$family)
  return(vcovObj)
}

#' @export
model.frame.glmmfit <- function(object, ...) {
    return(stats::model.frame(object$modObj))
}

#' Log-Likelihood for GLMM Fit
#'
#' This function extracts the log-likelihood of a fitted generalized linear mixed model (GLMM) object.
#'
#' @param object A fitted GLMM object of class \code{glmmfit}.
#' @param ... Additional arguments (currently not used).
#'
#' @return The log-likelihood value of the fitted model.
#'
#' @examples
#' \dontrun{
#'   fit <- glmmfit(data, family = "binomial")
#'   log_likelihood <- logLik(fit)
#' }
#'
#' @seealso \code{\link[stats]{logLik}} for the generic log-likelihood method.
#'
#' @export
logLik.glmmfit <- function(object, ...) {
    return(stats::logLik(object$modObj))
}


# beta <- c(3, 5)
# Sigma <- matrix(c(2,1,1,2),2)
# theta <- 0.39
# data <- batchGLMMData(beta,c(10,10),Sigma,100,rgen01(c(10,10)),"negative_binomial","log", theta=theta)
# fits <- glmmVpc::batchGLMMFit(Feature ~ X+(X|cluster), data, family="negative_binomial")
# test <- vpc.test(fits, Feature~1+(1|cluster))

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
  class(fits) <- "Glmmfits"
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


#' Extract Model Coefficients from a glmmfit Object
#'
#' This function extracts coefficients from a glmmfit object,
#' handling different model families (negative binomial, tweedie, gaussian).
#'
#' @param glmmfitObj An object of class glmmfit
#'
#' @return A named vector of coefficients. The exact components depend on the family:
#'   \itemize{
#'     \item For all families: Fixed effects coefficients (beta) and variance components (sigma_vec)
#'     \item For negative_binomial: Additional dispersion parameter (theta)
#'     \item For tweedie: Additional dispersion (phi) and power parameters
#'     \item For gaussian: Additional residual standard deviation (sigma_e)
#'   }
#'
#' @export
#' @method coef glmmfit
coef.glmmfit <- function(glmmfitObj) {
    family <- glmmfitObj$family
    beta <- glmmfitObj$beta
    names(beta) <- paste0("b", 1:length(beta))
    Sigma <- glmmfitObj$Sigma
    sigma_vec <- mat2vec(Sigma)
    switch(family,
           "negative_binomial" = {
              phi <- glmmfitObj$phi
              names(phi) <- "theta"
              c(beta,phi,sigma_vec)
           },
           "tweedie" = {
             phi <- glmmfitObj$phi
             power <- glmmfitObj$power
             names(phi) <- "phi"
             names(power) <- "power"
             c(beta,phi,sigma_vec,power)
           },
           "gaussian" = {
             sigma_e <- glmmfitObj$sigma_e
             names(sigma_e) <- "sigma_e"
             c(beta,sigma_e,sigma_vec)
           }
           )
}


#' Extract Model Coefficients from a Glmmfit Object
#'
#' This function extracts coefficients from a glmmfit object,
#' handling different model families (negative binomial, tweedie, gaussian).
#'
#' @param GlmmfitsObj An object of class glmmfit
#'
#' @return A named vector of coefficients. The exact components depend on the family:
#'   \itemize{
#'     \item For all families: Fixed effects coefficients (beta) and variance components (sigma_vec)
#'     \item For negative_binomial: Additional dispersion parameter (theta)
#'     \item For tweedie: Additional dispersion (phi) and power parameters
#'     \item For gaussian: Additional residual standard deviation (sigma_e)
#'   }
#'
#' @export
#' @method coef Glmmfits
coef.Glmmfits <- function(GlmmfitsObj) {
  result <- sapply(GlmmfitsObj, function(x) {
    return(stats::coef(x))
  }, simplify = TRUE)
  return(t(result))
}


