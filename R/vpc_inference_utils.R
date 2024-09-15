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
