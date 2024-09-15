#' parallelbatchGLMMData <- function(paramMat, X, ns, family, link = "identity", iter = 1, parallel = TRUE, num_cores = NULL, ...) {
#'   process_row <- function(i) {
#'     row <- as.numeric(paramMat[i, ])
#'     param_names <- names(paramMat)
#'     beta <- row[grep("b", param_names)]
#'     Sigma_vec <- row[grep("sig", param_names)]
#'     Sigma <- vec2mat(Sigma_vec)
#'
#'     # Dynamically capture additional parameters
#'     args <- list(...)
#'
#'     # Ensure family-specific parameters are passed correctly
#'     if (family == "negative_binomial" && is.null(args$theta)) {
#'       stop("'theta' must be specified for the Negative Binomial family.")
#'     }
#'     if (family == "tweedie" && (is.null(args$p) || is.null(args$phi))) {
#'       stop("'p' and 'phi' must be specified for the Tweedie family.")
#'     }
#'
#'     tryCatch({
#'       # Generate data using batchGLMMData with the family and link function
#'       gen_df <- batchGLMMData(num_feat = iter, X = X, beta = beta, Sigma = Sigma, ns = ns,
#'                               family = family, link = link, ...)
#'       return(gen_df[-(1:nrow(X)), ])  # Remove X rows (if X has multiple rows)
#'     }, error = function(e) {
#'       warning("Error in generating data for row ", i, ": ", conditionMessage(e))
#'       return(NULL)
#'     })
#'   }
#'
#'   if (parallel) {
#'     chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#'     if (nzchar(chk) && chk == "TRUE") {
#'       num_cores <- 2L
#'     } else {
#'       num_cores <- max(1, ifelse(is.null(num_cores), parallel::detectCores() - 1, num_cores))
#'     }
#'
#'     cl <- parallel::makeCluster(num_cores)
#'     on.exit(parallel::stopCluster(cl), add = TRUE)
#'
#'     # Export the necessary variables and functions to the worker processes
#'     parallel::clusterExport(cl, varlist = c("vec2mat", "batchGLMMData", "createGLMMData", "generateDataByFamily",
#'                                             "computeMuGLMM", "groups", "X", "mvn", "eta", "ns", "family", "link",
#'                                             "iter", "paramMat", "args"),
#'                             envir = environment())
#'
#'     dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row, cl = cl)
#'   } else {
#'     dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row)
#'   }
#'
#'   dataMatrix <- do.call(rbind, dataMatrix)
#'
#'   # Ensure row and column names are set correctly
#'   if (!is.null(dataMatrix) && nrow(dataMatrix) > 0) {
#'     rownames(dataMatrix) <- paste0("Feature", 1:nrow(dataMatrix))
#'     colnames(dataMatrix) <- groups(ns)
#'   }
#'
#'   # Return combined data (X + generated features)
#'   return(rbind(X, dataMatrix))
#' }
#'
