V_x <- function(mu, sigma_sq) {
  (exp(sigma_sq) - 1) * exp(2 * mu + sigma_sq)
}

Ex_p <- function(p, mu, sigma_sq) {
  exp(p * mu + (p^2 * sigma_sq) / 2)
}

vpc.nb <- function(mu, sigma_sq, phi){
  Vx <- V_x(mu, sigma_sq)
  Ex <- Ex_p(1, mu, sigma_sq)
  Ex2 <- Ex_p(2, mu, sigma_sq)
  Vx / (Vx + Ex + 1 / phi * Ex2)
}

vpc.tw  <- function(mu, sigma_sq, phi, p){
  Vx <- V_x(mu, sigma_sq)
  Exp <- Ex_p(p, mu, sigma_sq)
  Vx / (Vx + phi * Exp)
}

mu <- function(beta, x) {
  if(is.null(x)){
    beta
  } else {
    x <- c(1, x)
    sum(beta * x)
  }
}

sig <- function(sigma_vec, x) {

  z <- c(1, x)  # Ensure z has length matching Sigma dimensions
  if(length(sigma_vec) == 1) {
    sigma_vec
  }else if(length(sigma_vec) == 3 && length(x) == 1) {
    sig11 <- sigma_vec[1]
    sig12 <- sigma_vec[2]
    sig22 <- sigma_vec[3]
    sig11 + 2*sig12*x + sig22*x^2
  }

}

construct_beta_param_str <- function(beta) {
  n <- length(beta)
  beta_param_str <- paste0("b", 0:(n - 1), collapse = ", ")
  return(beta_param_str)
}

construct_sigma_param_str <- function(Sigma) {
  q <- ifelse(is.matrix(Sigma), nrow(Sigma), 1)
  Sigma_param_str <- paste0(unlist(sapply(1:q, function(k) {
    paste0('sig', k, k:q)
  })), collapse = ", ")
  return(Sigma_param_str)
}



##########    Schur's Complement  #######

A_block <- function(A, ind_lower) {
  n <- nrow(A)
  j <- ind_lower
  i <- setdiff(1:n,j)
  A11 <- A[i,i]
  A12 <- A[i,j]
  A22 <- A[j,j]

  return(list("A11"=A11,"A12"=A12,"A22"=A22))

}

schur_complement <- function(A, ind_lower) {
  partition <- A_block(A, ind_lower)

  A11 <- partition[["A11"]]
  A12 <- partition[["A12"]]; A21 <- t(A12)
  A22 <- partition[["A22"]]
  A11.inv <- solve(A11)
  return(A22 - A21 %*% A11.inv %*% A12)
}


#' Calculate mixing weight using Information Matrix Correlation
#'
#' This function computes mixing weight based on the correlation of parameters in the information matrix
#' of a model, using different methods including self-derived, Zhang, Julius, or all approaches.
#'
#' @param modObj A model object containing a Sigma matrix and for which vcov can be computed.
#' @param type A character string specifying the method to compute the p-value. Options are "self", "zhang", "julius", or "all". Default is "self".
#'
#' @return A numeric probability value or NULL if modObj$Sigma is not a matrix.
#'
#' @details
#' The function first checks if modObj$Sigma is a matrix. If not, it returns NULL.
#' It calculates the information matrix as the inverse of the variance-covariance matrix from vcov(modObj).
#'
#' - For type = "self", the correlation coefficient (rho) is derived directly from the off-diagonal and diagonal elements of the information matrix.
#' - For type = "zhang", "julius", or "all", the Schur complement of the information matrix subset is used to compute rho.
#'
#' The probability is computed as acos(rho) / (2 * pi), representing the angle of correlation in radians scaled to a probability.
#'
#' @examples
#' \dontrun{
#' # Assuming mod is a fitted model object
#' p_self <- pr11(mod, type = "self")
#' p_zhang <- pr11(mod, type = "zhang")
#' }
#' @export
pr11 <- function(modObj, type = c("self", "zhang", "julius", "all")) {
  type <- match.arg(type)
  if (!inherits(modObj$Sigma, "matrix")) return(NULL)
  inf.mat <- solve(stats::vcov(modObj))
  ind_lower <- grep("sig([0-9])\\1", colnames(inf.mat))
  #print(ind_lower)
  #print(modObj$modObj$sdr$pdHess)
  if (type == "self") {
    I_12 <- inf.mat[ind_lower[1], ind_lower[2]]
    I_11 <- inf.mat[ind_lower[1], ind_lower[1]]
    I_22 <- inf.mat[ind_lower[2], ind_lower[2]]
    rho <- I_12 /sqrt(I_11*I_22)
  } else if(type == "zhang" || type == "julius" || type == "all") {
    sch_comp <- schur_complement(inf.mat,ind_lower)
    rho <- sch_comp[1,2]/(sqrt(sch_comp[1,1]*sch_comp[2,2]))
  }
  p <- acos(rho)/(2*pi)
  if(is.na(p)) {
    #message("using ")
    p <- 0.25
  }
  return(p)
}



#' ############ q multivariate mix norm ########
#'
#' #' Compute Quantiles for a Multivariate Truncated Normal Distribution
#' #'
#' #' This function calculates the lower and upper quantiles for a multivariate
#' #' truncated normal distribution given the specified mean vector, covariance
#' #' matrix, and truncation limits.
#' #'
#' #' @param mu A numeric vector representing the mean of the distribution.
#' #' @param Sigma A covariance matrix for the multivariate normal distribution.
#' #' @param lower A numeric vector specifying the lower truncation limits.
#' #' @param upper A numeric vector specifying the upper truncation limits.
#' #' @param n_samples An integer specifying the number of samples to generate
#' #' (default is 10000).
#' #' @param alpha A numeric value representing the significance level for quantile
#' #' calculation (default is 0.05).
#' #'
#' #' @return A numeric vector containing the lower and upper quantiles corresponding
#' #' to the specified alpha.
#' #' @examples
#' #' mean_vector <- c(0, 0)
#' #' cov_matrix <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' #' lower_limits <- c(-Inf, -Inf)
#' #' upper_limits <- c(1, 1)
#' #' quantiles <- qmtmvnorm(mean_vector, cov_matrix, lower_limits, upper_limits)
#' #' print(quantiles)
#' #' @export
#' qmtmvnorm <- function(mu, Sigma, lower, upper, n_samples = 10000, alpha = 0.05) {
#'   samples <- rtmvnorm(n_samples, mean = mu, sigma = Sigma,
#'                       lower = lower, upper = upper)
#'
#'   quantiles <- quantile(samples, probs = c(alpha / 2, 1 - alpha / 2))
#'   return(quantiles)
#' }
