#' singleGLMMFit <- function(formula, data, family, cov_values, nulltest=NULL) {
#'   glmmTMBfamily <- glmmTMBfamily(family)
#'   modObj <- tryCatch({
#'     mod <- glmmTMB::glmmTMB(formula=formula, data=data, family=glmmTMBfamily)
#'     mod
#'   }, error = function(e) {
#'     message("Error fitting model: ", e$message)
#'     return(NULL)
#'   })
#'
#'   if (is.null(modObj)) {
#'     return(NULL)  # Return NULL if the model fitting failed
#'   }
#'
#'   beta <- unname(glmmTMB::fixef(modObj)$cond)
#'   names(beta) <- paste0("b", 0:(length(beta)-1))
#'
#'   Sigma <- unname(glmmTMB::VarCorr(modObj)$cond)[[1]]
#'
#'   phi <- glmmTMB::sigma(modObj)
#'   names(phi) <- "phi"
#'
#'   power <- unname(glmmTMB::family_params(modObj))
#'
#'
#'   if (!identical(power, numeric(0))) names(power) <- "power"
#'
#'   vcov_ <- tryCatch({
#'     vcov.fit(modObj)
#'   }, error = function(e) {
#'     message("Error extracting variance-covariance matrix: ", e$message)
#'     return(NULL)
#'   })
#'
#'   if (is.null(vcov_)) {
#'     return(NULL)  # Return NULL if variance-covariance extraction failed
#'   }
#'
#'   n <- modObj$modelInfo$nobs
#'
#'   vpcfunc_name <- paste0("vpc.", family)
#'   vpcfunc <- get(vpcfunc_name, mode="function", inherits=TRUE)
#'
#'   vpc <- sapply(cov_values, function(val) {
#'     v <- vpcfunc(beta=beta, Sigma=Sigma, power=power, phi=phi, x=val)
#'     names(v) <- paste0("vpc", val)
#'     return(v)
#'   })
#'
#'   res <- list("beta"=beta, "Sigma"=Sigma, "phi"=phi, "power"=power,
#'               "family"=family, "vcov" = vcov_, "n" = n, "vpc"=vpc,
#'               "convergence"=modObj$fit$convergence,"message"=modObj$fit$message,
#'               "pdHess"=modObj$sdr$pdHess)
#'   class(res) <- "glmmfit"
#'   return(res)
#' }
#'
#'
#' #' Fit Multiple Generalized Linear Mixed Models (GLMMs)
#' #'
#' #' This function fits multiple Generalized Linear Mixed Models (GLMMs) to the rows of a data matrix.
#' #' Each row of the matrix is treated as a separate dataset, and the function returns a list of fitted models.
#' #'
#' #' @param formula A formula specifying the fixed and random effects in the model.
#' #' @param dataMat A matrix where each row represents a separate dataset to fit the model to.
#' #' @param X A numeric vector or matrix of covariates for the fixed effects.
#' #' @param group A factor or numeric vector indicating the grouping structure for random effects.
#' #' @param family A character string specifying the family of distributions to use (e.g., "nb" for negative binomial).
#' #' @param cov_values A vector of covariate values at which to evaluate the Variance Partition Coefficient (VPC).
#' #' @param nulltest An optional argument for specifying a null hypothesis test (default is NULL).
#' #' @param parallel Logical, whether to run the fitting in parallel. Default is TRUE.
#' #' @param num_cores Integer, the number of cores to use for parallel processing. Defaults to the number of available cores minus one.
#' #'
#' #' @return A list containing the following components:
#' #' \item{GlmmObj}{A list of GLMM objects, one for each row of the input data matrix.}
#' #' \item{num_feat}{The number of datasets (rows of the input matrix) that were fit.}
#' #' \item{num_covariates}{The number of covariate values used for VPC calculations.}
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' ns <- c(10, 20)
#' #' X <- rgen01(ns)
#' #' beta <- c(0.5, 2)
#' #' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
#' #' family <- "nb"
#' #' cov_values <- c(0.5, 1.5)
#' #' formula <- y ~ x + (1 | group)
#' #' datafrmMat <- batchGLMMData(5, X, beta, Sigma, ns, family = "nb", theta = 1.5)
#' #' ys <- datafrmMat[-1,]; X <- datafrmMat[1,]; group <- colnames(datafrmMat)
#' #' fit <- batchGLMMFit(formula, ys, X, group, family, cov_values)
#' batchGLMMFit <- function(formula, dataMat, X, group, family, cov_values, nulltest=NULL, parallel=TRUE, num_cores=NULL) {
#'
#'   if (parallel) {
#'
#'     chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
#'     if (nzchar(chk) && chk == "TRUE") {
#'       num_cores <- 2L
#'     } else {
#'       num_cores <-ifelse(!is.null(num_cores), num_cores,parallel::detectCores() - 1)
#'     }
#'
#'     cl <- parallel::makeCluster(num_cores)
#'     on.exit(parallel::stopCluster(cl), add = TRUE)
#'
#'     parallel::clusterExport(cl, varlist=c("singleGLMMFit", "formula", "X", "group", "family", "cov_values", "nulltest"), envir=environment())
#'
#'     fitobj <- pbapply::pblapply(1:nrow(dataMat), function(i){
#'       data <- data.frame(x=X, y=dataMat[i,], group=group)
#'       fit <- singleGLMMFit(formula = formula, data=data, family=family, cov_values=cov_values, nulltest = nulltest)
#'       return(fit)
#'     }, cl = cl)
#'
#'   } else {
#'     fitobj <- pbapply::pblapply(1:nrow(dataMat), function(i){
#'       data <- data.frame(x=X, y=dataMat[i,], group=group)
#'       fit <- singleGLMMFit(formula = formula, data=data, family=family, cov_values=cov_values, nulltest = nulltest)
#'       return(fit)
#'     })
#'   }
#'
#'   #fitobj <- fitobj[!sapply(fitobj, is.null)]
#'   names(fitobj) <- paste0("Feature", 1:nrow(dataMat))
#'   fitobj <- list("GlmmObj"=fitobj, "num_feat"=nrow(dataMat), "num_covariates"=length(cov_values))
#'   class(fitobj) <- "Glmmfit"
#'   return(fitobj)
#' }
#'
#'
#' #' Extract Coefficients from a GLMM Fit Object
#' #'
#' #' This function extracts key parameters such as `beta`, `Sigma`, `phi`, `power`, and `vpc`
#' #' from a Generalized Linear Mixed Model (GLMM) fit object and returns them as a single vector.
#' #'
#' #' @param fitObj An object containing the fit results from a GLMM, with components `beta`, `Sigma`, `phi`, `power`, and `vpc`.
#' #'
#' #' @return A numeric vector containing the coefficients `beta`, `phi`, the vectorized version of the covariance matrix `Sigma`, `power`, and `vpc`.
#' #' @export
#' #'
#' #' @examples
#' coef.glmmfit <- function(object, ...) {
#'
#'   beta <- object$beta
#'   Sigma <- object$Sigma
#'   phi <- object$phi
#'   power <- object$power
#'   vpc <- object$vpc
#'   c(beta, phi, mat2vec(Sigma), power, vpc)
#' }
#'
#' #' Extract Coefficients from Glmmfit Object
#' #'
#' #' This function extracts coefficients from a batch of GLMM fit objects stored within a `Glmmfit` object.
#' #'
#' #' @param object A `Glmmfit` object that contains multiple GLMM fit results.
#' #' @param ... Additional arguments (currently not used).
#' #'
#' #' @return A data frame where each row represents the coefficients from one of the fitted GLMM models.
#' #' @export
#' #' @importFrom stats coef
#' #'
#' #' @examples
#' #' ns <- c(10, 20)
#' #' X <- rgen01(ns)
#' #' beta <- c(0.5, 2)
#' #' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
#' #' family <- "nb"
#' #' cov_values <- c(0.5, 1.5)
#' #' formula <- y ~ x + (1 | group)
#' #' datafrmMat <- batchGLMMData(5, X, beta, Sigma, ns, family = "nb", theta = 1.5)
#' #' ys <- datafrmMat[-1,]; X <- datafrmMat[1,]; group <- colnames(datafrmMat)
#' #' fit <- batchGLMMFit(formula, ys, X, group, family, cov_values)
#' #' coefs <- coef(fit)
#' coef.Glmmfit <- function(object, ...) {
#'   GlmmObj <- object$GlmmObj
#'   num_feat <- object$num_feat
#'
#'   # Extract coefficients from each GLMM fit
#'   coef_lists <- lapply(1:num_feat, function(x) {
#'     text=paste0("GlmmObj$Feature",x)
#'     coef(eval(parse(text=text)))
#'   })
#'
#'   # Combine the coefficients into a data frame
#'   coefs <- data.frame(do.call(rbind, coef_lists))
#'   rownames(coefs) <- paste0("Feature", which(!sapply(coef_lists, is.null)))
#'
#'   return(coefs)
#' }
#'
#'
#' #' Variance Partition Coefficient for Negative Binomial Model
#' #'
#' #' This function calculates the Variance Partition Coefficient (VPC)
#' #' for a negative binomial model, given the parameters.
#' #'
#' #' @param beta A vector of regression coefficients.
#' #' @param Sigma A covariance matrix for the random effects.
#' #' @param phi The overdispersion parameter (can be overridden by `theta`).
#' #' @param x The covariate values.
#' #' @param ... Optional additional arguments, such as `theta` to override `phi`.
#' #'
#' #' @return The VPC for the negative binomial model.
#' #' @export
#' #'
#' #' @examples
#' #' beta <- c(0.5, 1.2)
#' #' Sigma <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2)
#' #' phi <- 1.5
#' #' x <- 0.7
#' #' vpc.nb(beta, Sigma, phi, x)
#' vpc.nb <- function(beta, Sigma, phi, x,...) {
#'   args <- list(...)
#'   if (!is.null(args$theta)) phi <- args$theta
#'   X <- Z <- c(1, x)
#'   if (length(beta) != length(X)) {
#'     stop("Length of beta must match the number of elements in X.")
#'   }
#'
#'   # if (dim(Sigma)[1] != length(Z)) {
#'   #   stop("row of Sigma matrix dimensions must match the length of Z.")
#'   # }
#'   mu <- as.numeric(X %*% beta)
#'   sig <- as.numeric(Z %*% Sigma %*% Z)
#'   Ex <- logNormM(mu,sig,1)
#'   Ex.2 <- logNormM(mu,sig,2)
#'   Vx <- Ex.2 - Ex^2
#'   vpcnb <- Vx/(Vx+Ex+(1/phi)*Ex.2)
#'   vpcnb
#' }
#'
#' #' Variance Partition Coefficient for Tweedie (Compound Poisson)
#' #'
#' #' This function calculates the Variance Partition Coefficient (VPC)
#' #' for a transformed Weibull model with a power transformation.
#' #'
#' #' @param beta A vector of regression coefficients.
#' #' @param Sigma A covariance matrix for the random effects.
#' #' @param phi The overdispersion parameter.
#' #' @param power The power parameter.
#' #' @param x The covariate values.
#' #' @param ... Optional additional arguments.
#' #'
#' #' @return The VPC for the transformed Weibull model.
#' #' @export
#' #'
#' #' @examples
#' #' beta <- c(0.5, 1.2)
#' #' Sigma <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2)
#' #' phi <- 1.5
#' #' power <- 2
#' #' x <- 0.7
#' #' vpc.tw(beta, Sigma, phi, power, x)
#' vpc.tw <- function(beta, Sigma, phi, power, x,...) {
#'   X <- Z <- c(1, x)
#'   if (length(beta) != length(X)) {
#'     stop("Length of beta must match the number of elements in X.")
#'   }
#'
#'   # if (dim(Sigma)[1] != length(Z)) {
#'   #   stop("Sigma matrix dimensions must match the length of Z.")
#'   # }
#'   mu <- X %*% beta
#'   sig <- Z %*% Sigma %*% Z
#'   Ex <- logNormM(mu,sig,1)
#'   Ex2 <- logNormM(mu,sig,2)
#'   Ex.p <- logNormM(mu,sig,power)
#'   Vx <- Ex2 - Ex^2
#'   vpctw <- Vx/(Vx+phi*Ex.p)
#'   vpctw
#' }
#'
#' #' Variance Partition Coefficient for Gaussian Model
#' #'
#' #' This function calculates the Variance Partition Coefficient (VPC)
#' #' for a Gaussian model, adjusting for random effects and overdispersion.
#' #'
#' #' @param Sigma A covariance matrix for the random effects.
#' #' @param phi The overdispersion parameter.
#' #' @param x The covariate values.
#' #' @param ... Optional additional arguments.
#' #'
#' #' @return The VPC for the Gaussian model.
#' #' @export
#' #'
#' #' @examples
#' #' Sigma <- matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2)
#' #' phi <- 1.5
#' #' x <- 0.7
#' #' vpc.ga(Sigma, phi, x)
#' vpc.ga <- function(Sigma, phi, x, ...) {
#'   Z <- c(1, x)
#'   # if (dim(Sigma)[1] != length(Z)) {
#'   #   stop("Sigma matrix dimensions must match the length of Z.")
#'   # }
#'   sig <- Z %*% Sigma %*% Z
#'   sig/(sig+phi^2)
#' }
