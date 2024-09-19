#' Generate GLMM Data for Multiple Parameter Sets in Parallel
#'
#' This function generates generalized linear mixed model (GLMM) data for multiple parameter sets.
#' It allows parallel processing for faster computation by distributing the tasks across multiple cores.
#'
#' @param paramMat A matrix where each row represents a set of parameters for generating GLMM data.
#'        Columns should include fixed effect coefficients (beta) and random effect covariances (Sigma).
#' @param X A matrix of predictors (design matrix) for the fixed effects.
#' @param ns A vector specifying the number of observations per group (cluster sizes).
#' @param family A string specifying the distribution family. Supported families include
#'               "gaussian", "poisson", "negative_binomial", "tweedie", and "gamma".
#' @param link A string specifying the link function (default is "identity").
#' @param iter An integer specifying the number of features to generate per parameter set (default is 1).
#' @param parallel A logical value indicating whether to run the computation in parallel (default is TRUE).
#' @param num_cores The number of cores to use for parallel processing. If NULL, the function will use all
#'                  available cores minus one (default is NULL).
#'
#' @export
#'
#' @return A matrix where rows represent the predictors from `X` followed by the generated feature data.
#'         If parallel processing is used, the data is generated in parallel across multiple cores.
batchGLMMDataUsingMatrix <- function(paramMat, ns, X=NULL,
                                  family="gaussian",
                                  link = "identity",
                                  iter = 1, parallel = TRUE,
                                  num_cores = NULL) {

  if (parallel && is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 1  # Reserve one core for system tasks
  }

  # Helper function to process each row of paramMat
  process_row <- function(i) {
    row <- paramMat[i, ]
    family_args <- validate_family_params(family,row)
    row <- as.numeric(row)
    param_names <- names(paramMat)
    beta <- row[grep("b", param_names)]
    Sigma_vec <- row[grep("sig", param_names)]
    Sigma <- vec2mat(Sigma_vec)
    tryCatch({
      gen_df <- do.call(batchGLMMData, c(list(
        beta = beta,
        ns = ns,
        Sigma = Sigma,
        num = iter,
        X = X,
        family = family,
        link = link),
        family_args ))
      gen_df <- get_y(gen_df)
      return(gen_df)
    }, error = function(e) {
      warning("Error in generating data for row ", i, ": ", conditionMessage(e))
      return(NULL)
    })

  }

  print(paste("Generating Data from ", family))
  if (parallel) {
    dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row, cl = num_cores)
  } else {
    dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row)
  }

  dataMatrix <- do.call(rbind, dataMatrix)

  #Ensure proper row and column names
  if (!is.null(dataMatrix) && nrow(dataMatrix) > 0) {
    rownames(dataMatrix) <- paste0("Feature", 1:nrow(dataMatrix))
    colnames(dataMatrix) <- cluster_assignment(ns)
  }

  if (!is.null(X)) {
    Xt <- t(X)
    if (nrow(Xt) == 1) {
      rownames(Xt) <- "X"
    } else {
      rownames(Xt) <- paste0("X", 1:nrow(Xt))
    }
  }

  fullData <- rbind(Xt, dataMatrix)
  class(fullData) <- c("batchglmmDataMatrix", "matrix", "array")
  attr(fullData, "num_feat") <- nrow(dataMatrix)
  attr(fullData, "num_covariate") <- ifelse(is.null(X), 0, nrow(Xt))

  return(fullData)
}




# parallelbatchGLMMData <- function(paramMat, X, ns, family, link = "identity",
#                                   iter = 1, parallel = TRUE,
#                                   num_cores = NULL, ...) {
#
#   # Validate family parameters once before starting the batch generation
#   args <- list(...)
#   validate_family_params(family, args)  # Call the validation function
#
#   process_row <- function(i) {
#     row <- as.numeric(paramMat[i, ])
#     param_names <- names(paramMat)
#
#     # Extract beta and Sigma from paramMat
#     beta <- row[grep("b", param_names)]
#     Sigma_vec <- row[grep("sig", param_names)]
#     Sigma <- vec2mat(Sigma_vec)
#     # args <- list(...)
#
#     tryCatch({
#       gen_df <- batchGLMMData(num_feat = iter, X = X, beta = beta,
#                               Sigma = Sigma, ns = ns,
#                               family = family, link = link, ...)
#
#       return(gen_df[-(1:nrow(X)), ])   # Remove X rows
#     }, error = function(e) {
#       warning("Error in generating data for row ", i, ": ", conditionMessage(e))
#       return(NULL)
#     })
#   }
#
#   # Determine the number of cores for parallel processing
#   if (parallel) {
#     chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#     if (nzchar(chk) && chk == "TRUE") {
#       num_cores <- 2L
#     } else {
#       num_cores <- max(1, ifelse(is.null(num_cores),
#                                  parallel::detectCores() - 1, num_cores))
#     }
#
#     cl <- parallel::makeCluster(num_cores)
#     on.exit(parallel::stopCluster(cl), add = TRUE)
#
#     # Export the necessary variables and functions to the worker processes
#     parallel::clusterExport(cl, varlist = c("vec2mat", "batchGLMMData",
#                                             "singleGLMMData",
#                                             "generateDataByFamily",
#                                             "computeMuGLMM",
#                                             "generateRandomEffectsCoefficients",
#                                             "cluster_assignment",
#                                             "X", "eta", "ns",
#                                             "family", "link", "iter",
#                                             "paramMat", "args"),
#                             envir = environment())
#
#     dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row, cl = cl)
#   } else {
#     dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row)
#   }
#
#   # Combine the results into a single matrix
#   dataMatrix <- do.call(rbind, dataMatrix)
#
#   # Ensure row and column names are set correctly
#   if (!is.null(dataMatrix) && nrow(dataMatrix) > 0) {
#     rownames(dataMatrix) <- paste0("Feature", 1:nrow(dataMatrix))
#     colnames(dataMatrix) <- groups(ns)
#   }
#   return(rbind(X, dataMatrix))
# }
#
# validate_family_params <- function(family, args) {
#   # Define required parameters for each family
#   required_params <- list(
#     "negative_binomial" = c("theta"),
#     "tweedie" = c("p", "phi"),
#     "poisson" = c(),  # No additional parameters
#     "gaussian" = c()  # No additional parameters
#     # Add more families and their required parameters here
#   )
#
#   # Check if the family is recognized
#   if (!family %in% names(required_params)) {
#     stop(paste("Unsupported family:", family))
#   }
#
#   # Get the required parameters for the given family
#   required <- required_params[[family]]
#
#   # Check if the required parameters are present in `args`
#   missing_params <- setdiff(required, names(args))
#   if (length(missing_params) > 0) {
#     stop(paste("Missing required parameters for", family, "family:", paste(missing_params, collapse = ", ")))
#   }
#
#   # If all parameters are present, return TRUE
#   return(TRUE)
# }



