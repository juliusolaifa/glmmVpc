#' Calculate VPC for a specific family
#'
#' This function computes the Variance Partition Coefficient (VPC) for a given
#' model family (e.g., Negative Binomial, Tweedie, Gaussian, comPoisson).
#' The necessary parameters like `beta`, `Sigma`, and family-specific parameters
#' (e.g., `phi`, `power`, `nu`) are passed via the `args` list.
#'
#' @param family A string representing the family of the model (e.g.,
#' "negative_binomial", "tweedie", "gaussian", "compoisson").
#' @param args A list containing parameters for the VPC calculation.
#' Must include `beta` (regression coefficients),
#'        `Sigma` (covariance matrix for the random effects), and
#'        family-specific parameters:
#'        \itemize{
#'          \item `phi` for Negative Binomial and Tweedie families.
#'          \item `power` for Tweedie family.
#'          \item `nu` for comPoisson family.
#'        }
#' @param x The covariate values used to calculate `mu` and `sig`.
#'
#' @return The calculated Variance Partition Coefficient (VPC) for
#' the specified family.
#' @export
#'
#' @examples
#' # Example args list for Tweedie family
#' args <- list(
#'   beta = c(0.5, 1.2),
#'   Sigma = matrix(c(0.1, 0.05, 0.05, 0.2), nrow = 2),
#'   phi = 1.5,
#'   power = 1.6
#' )
#'
#' # Calculate VPC for Tweedie family
#' vpc_value <- calculate_vpc_for_family(family = "tweedie", args = args, x = 0.5)
#' print(vpc_value)
calculate_vpc_for_family <- function(family, args, x) {

  beta <- args$beta
  Sigma <- args$Sigma
  if (is.null(beta)) stop("beta parameter must be provided for vpc.")
  if (is.null(Sigma)) stop("Sigma parameter must be provided for vpc.")
  result <- compute_mu_sig(beta=beta, Sigma=Sigma, x=x)
  mu <- result$mu
  sig <- result$sig

  vpc_value <- switch(family,
                      "negative_binomial" = {
                        phi <- args$phi
                        if (is.null(phi)) stop("phi parameter must be
                                                 provided for negative binomial
                                               family.")
                        moments <- compute_lnmoments(mu, sig)
                        Vx <- moments$Ex2 - moments$Ex^2
                        Vx / (Vx + moments$Ex + (1 / phi) * moments$Ex2)
                      },
                      "tweedie" = {
                        power <- args$power
                        phi <- args$phi
                        if (is.null(phi)) stop("phi parameter must be
                                                 provided for tweedie
                                               family.")
                        if (is.null(power)) stop("Power parameter must be
                                                 provided for Tweedie family.")
                        moments <- compute_lnmoments(mu, sig, power)
                        Vx <- moments$Ex2 - moments$Ex^2
                        Vx / (Vx + phi * moments$Ex_p)
                      },
                      "gaussian" = {
                        if(!is.null(x)) {
                          z <- c(1,x)
                          sig <- as.numeric(z %*% Sigma %*% z)
                        } else{
                          sig <- as.numeric(Sigma)
                        }
                        sigma_e <- args$sigma_e
                        sig / (sig + sigma_e^2)
                      },
                      "compoisson" = {
                        nu <- args$nu
                        if (is.null(nu)) stop("nu parameter must be provided for
                                              comPoisson family.")
                        moments <- compute_lnmoments(mu, sig)
                        Vx <- moments$Ex2 - moments$Ex^2
                        Vx / (Vx + moments$Ex + (1 / nu) * moments$Ex2)
                      },
                      stop("Unsupported family for VPC calculation")
  )

  return(vpc_value)
}


#' Variance Partition Coefficient (VPC) Calculation for GLMM Models
#'
#' This function calculates the Variance Partition Coefficient (VPC)
#' for a fitted GLMM model, regardless of the family.
#'
#' @param model_fit A fitted GLMM model object from `singleGLMMFit`.
#' @param x The covariate values for which to calculate VPC.
#' @param ... Additional arguments for family-specific parameters
#' (e.g., `power` for Tweedie, `nu` for comPoisson).
#'
#' @return The Variance Partition Coefficient (VPC) for the specified model.
#' @export
vpc <- function(model_fit, x=NULL, ...) {

  calculate_vpc_single <- function(single_fit) {
    family <- single_fit$family
    vpc_value <- calculate_vpc_for_family(family=family, args=single_fit, x=x)

    #vcov_matrix <- stats::vcov(single_fit)
    #n <- stats::nobs(single_fit)

    result <- list(vpc_value,
              modObj = single_fit
      # beta = single_fit$beta,
      # Sigma = single_fit$Sigma,
      # phi = single_fit$phi,
      # p = single_fit$p,
      # family = family,
      ,x = x
      # vcov = vcov_matrix,
      # n = n
      )

    class(result) <- "vpcObj"
    return(result)
  }

  if(inherits(model_fit, "glmmfit")) {
    return(calculate_vpc_single(model_fit))
  }

  if (is.list(model_fit)) {
    if (all(sapply(model_fit, inherits, "glmmfit"))) {
      result <- (lapply(model_fit, calculate_vpc_single))
      class(result) <- "VpcObj"
      return(result)
    }
  }

  stop("Invalid input: model_fit should be of class 'glmmfit' or a list of 'glmmfit' models.")

}

#' Calculate Variance Partition Coefficient (VPC) for Different Families
#'
#' This function takes a family type, a parameter matrix, and a covariate
#' (or covariates),
#' and calculates the VPC for each row of the parameter matrix. The function
#' dynamically
#' handles family-specific parameters such as \code{phi} and \code{power} for Tweedie
#' and Negative Binomial families.
#'
#' @param family A character string specifying the family for the model.
#' Can be one of
#'   \code{"negative_binomial"} or \code{"tweedie"}.
#' @param paramMat A matrix or data frame containing the parameters for each
#' row. Columns
#'   for \code{b0}, \code{b1}, \code{sig11}, \code{sig12}, \code{sig22},
#'   and family-specific
#'   parameters (such as \code{theta} for Negative Binomial or \code{phi},
#'   \code{power} for Tweedie)
#'   are required.
#' @param x Covariates or other necessary inputs for calculating the VPC.
#'
#' @return A list where each element corresponds to the VPC results for each row
#' of the parameter matrix.
#' @export
vpc_from_paramMat <- function(family, paramMat, x) {
  vpc_results <- vector("list", nrow(paramMat))
  beta_cols <- grep("^b", colnames(paramMat))
  sigma_cols <- grep("^sig", colnames(paramMat))

  for (i in seq_len(nrow(paramMat))) {
    row <- as.list(paramMat[i, ])
    beta <- unlist(row[beta_cols])
    sigma_vec <- unlist(row[sigma_cols])
    Sigma <- vec2mat(sigma_vec)

    args <- list(
      beta = beta,
      Sigma = Sigma
    )
    args <- switch(
      family,
      "negative_binomial" = c(args, list(phi = row$theta)),
      "tweedie" = c(args, list(phi = row$phi, power = row$power)),
      args
    )

    vpc_results[[i]] <- calculate_vpc_for_family(family, args, x)
  }
  return(vpc_results)
}


