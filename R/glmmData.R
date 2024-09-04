#' Generate Data from a Generalized Linear Mixed Model (GLMM)
#'
#' This function generates simulated data from a Generalized Linear Mixed Model
#' (GLMM) using the specified fixed effects, random effects, and family of
#' distributions. It supports negative binomial and Tweedie distributions.
#'
#' @param X A numeric vector of fixed effect covariates.
#' @param beta A numeric vector of fixed effect coefficients. The length of
#' `beta` should match the length of `X`.
#' @param Sigma A numeric matrix representing the covariance matrix for the
#' random effects.
#' @param ns A numeric vector specifying the number of observations per group.
#' @param family A character string specifying the family of distributions to
#' use. Supported values are `"nb"` for negative binomial and `"tw"` for
#' Tweedie distribution.
#' @param ... Additional arguments passed to the family-specific random number
#' generation functions. For `"nb"` (negative binomial), use `theta`. For `"tw"`
#' (Tweedie), use `phi` and `power`.
#'
#' @return A numeric vector of generated data based on the specified GLMM.
#'
#' @export
#'
#' @examples
#' ns <- c(10, 20)
#' X <- rgen01(ns)
#' beta <- c(0.5, 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
#' data_nb <- createGLMMData(X, beta, Sigma, ns, family = "nb", theta = 1.5)
#' data_tw <- createGLMMData(X, beta, Sigma, ns, family = "tw", phi = 1.2, power = 1.5)
createGLMMData <- function(X, beta, Sigma, ns, family, ...) {
  ng <- length(ns)
  u <- mvn(ng, Sigma)
  X <- stats::model.matrix(~X)
  z0 <- stats::model.matrix(~ 0 + groups(ns))
  Z <- do.call(cbind, lapply(1:ncol(u), function(i) z0 * X[, i]))
  eta <-  eta(X, beta, Z, c(u))
  mu = split(x=exp(eta), f = rep(1:length(ns), ns))
  args <- list(...)
  data <- switch(family,
                 nb=mapply(MASS::rnegbin, n=ns, mu=mu, theta=args$theta),
                 tw=mapply(tweedie::rtweedie, n=ns, mu=mu, phi=args$phi, power=args$power),
                 stop(family, ": not implemented")
  )
  data <- if (methods::is(data, "list")) do.call(c, data) else c(data)
}


#' Generate Multiple Datasets from a Generalized Linear Mixed Model (GLMM)
#'
#' This function generates multiple datasets from a Generalized Linear Mixed
#' Model (GLMM) using the specified fixed effects, random effects, and family
#' of distributions. The generated datasets are combined into a single matrix.
#'
#' @param num_feat An integer specifying the number of datasets (features)
#' to generate.
#' @param X A numeric vector of fixed effect covariates.
#' @param beta A numeric vector of fixed effect coefficients. The length of
#' `beta` should match the length of `X`.
#' @param Sigma A numeric matrix representing the covariance matrix for the
#' random effects.
#' @param ns A numeric vector specifying the number of observations per group.
#' @param family A character string specifying the family of distributions to
#' use. Supported values are `"nb"` for negative binomial and `"tw"` for
#' Tweedie distribution.
#' @param ... Additional arguments passed to the family-specific random number
#' generation functions. For `"nb"` (negative binomial), use `theta`. For `"tw"`
#' (Tweedie), use `phi` and `power`.
#'
#' @return A numeric matrix where each row represents a generated dataset and
#' each column corresponds to a group specified in `ns`. The first rows of the
#' matrix correspond to the original `X` vector.
#'
#' @export
#'
#' @examples
#' ns <- c(10, 20)
#' X <- rgen01(ns)
#' beta <- c(0.5, 2)
#' Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
#' batchGLMMData(5, X, beta, Sigma, ns, family = "nb", theta = 1.5)
#' batchGLMMData(5, X, beta, Sigma, ns, family = "tw", phi = 1.2, power = 1.5)
batchGLMMData <- function(num_feat, X, beta, Sigma, ns, family, ...) {
  dataMatrix <- lapply(1:num_feat, function(x){
    createGLMMData(X, beta, Sigma, ns, family, ...)
  })
  dataMatrix <- do.call(rbind, dataMatrix)
  rownames(dataMatrix) <- paste("Feature ", 1:nrow(dataMatrix))
  colnames(dataMatrix) <- groups(ns)
  return(rbind(X, dataMatrix))
}



#' Generate GLMM Data from a Parameter Matrix
#'
#' This function generates datasets from a Generalized Linear Mixed Model (GLMM)
#' using a matrix of parameters. The function supports parallel processing
#' to speed up the generation process, especially with large parameter matrices.
#'
#' @param paramMat A matrix where each row contains parameters for generating a GLMM dataset.
#' @param X A numeric vector of fixed effect covariates.
#' @param ns A numeric vector specifying the number of observations per group.
#' @param family A character string specifying the family of distributions to use.
#' Supported values are `"nb"` for negative binomial and `"tw"` for Tweedie distribution.
#' @param iter An integer specifying the number of datasets (features) to generate per parameter row.
#' @param parallel A logical value indicating whether to use parallel processing. Default is `TRUE`.
#' @param num_cores integer indicating number of cores to use fo fitting
#'
#' @return A numeric matrix where each row represents a generated dataset, and each column corresponds to a group specified in `ns`.
#' @export
#'
#' @examples
#' ns <- c(10, 20)
#' X <- rgen01(ns)
#' family <- "nb"
#' b0 <- 3
#' b1 <- 7
#' phi <- 0.8
#' sig11 <- 2
#' sig12 <- 1
#' sig22 <- 2
#' paramMat <- expand.grid(b0=b0,b1=b1,phi=phi,sig11=sig11,sig12=sig12,sig22=sig22)
#' parallelbatchGLMMData(paramMat, X, ns, family, iter=3, parallel=TRUE)
parallelbatchGLMMData <- function(paramMat, X, ns, family, iter=1, parallel=TRUE, num_cores=NULL) {

  process_row <- function(i) {
    row <- as.numeric(paramMat[i, ])
    param_names <- names(paramMat)
    beta <- row[grep("b", param_names)]
    Sigma_vec <- row[grep("sig", param_names)]
    Sigma <- vec2mat(Sigma_vec)
    phi <- row[grep("phi", param_names)]
    power <- row[grep("power", param_names)]

    gen_df <- batchGLMMData(num_feat=iter, X=X, beta=beta, Sigma=Sigma, ns=ns,
                       family=family, theta=phi, power=power)
    return(gen_df[-1,])
  }

  if (parallel) {

    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      num_cores <- 2L
    } else {
      num_cores <- parallel::detectCores() - 1
    }

    cl <- parallel::makeCluster(num_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(cl, varlist=c("vec2mat", "batchGLMMData", "createGLMMData", "groups", "X",
                                "mvn", "eta", "ns", "family", "iter", "paramMat"),
                  envir=environment())

    dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row, cl = cl)
  } else {
    dataMatrix <- pbapply::pblapply(1:nrow(paramMat), process_row)
  }

  dataMatrix <- do.call(rbind, dataMatrix)
  total_rows <- nrow(paramMat) * iter
  rownames(dataMatrix) <- paste0("Feature", 1:total_rows)
  colnames(dataMatrix) <- groups(ns)
  return(rbind(X, dataMatrix))
}
