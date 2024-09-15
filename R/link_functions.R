#' Compute Mu from Eta using Link Function
#'
#' This function applies the specified link function to the linear predictor
#' `eta` to compute the mean `mu`. Supported link functions include identity,
#' log, logit, and inverse.
#'
#' @param eta The linear predictor. Must be a numeric vector.
#' @param link The link function to apply. Supported links are "identity",
#' "log", "logit", and "inverse".
#'
#' @return The mean `mu` after applying the link function.
#' @export
#'
#' @examples
#' eta <- c(0.5, 1.2, -0.7)
#' computeMuFromEta(eta, "logit")
#'
computeMuFromEta <- function(eta, link = "identity") {
  if (!is.numeric(eta)) {
    stop("Input 'eta' must be a numeric vector.")
  }
  if (any(is.na(eta))) {
    warning("NA values detected in 'eta'. Results will contain NAs.")
  }
  switch(link,
         "identity" = eta,
         "log" = exp(eta),
         "logit" = 1 / (1 + exp(-eta)),
         "inverse" = {
           if (any(eta == 0)) {
             stop("Zero values in 'eta' are not allowed with the 'inverse' link function.")
           }
           1 / eta
         },
         stop(paste("Unsupported link function:", link))
  )
}
